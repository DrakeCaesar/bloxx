#include <algorithm>
#include <ankerl/unordered_dense.h>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include "config.h"
#include "display.h"
#include "grid_config.h"
#include "search.h"
#include "search_bfs.h"
#include "search_dfs.h"
#include "search_impl.h"
#include "state_utils.h"
#include "transition.h"

using namespace std;

// When non-null and open, baseline build path is written here instead of cout
ofstream *g_baselineFile = nullptr;

// Number of lines last drawn for "New best" block (for in-place redraw)
static int s_newBestBlockLines = 0;

// Globals (declared in search_impl.h)
int g_numPos;
int g_neighborShifts[49][4];
int g_neighborCount[49];
int g_maxType[49];
StateType g_baseState;
const uint32_t g_requiredMask[6] = {0, 0, 2, 6, 14, 30};
int g_numSymm;
int g_symmPerm[8][49];

static void buildSymmetries() {
  int R = GRID.rows, C = GRID.cols;
  // All 8 elements of the dihedral group D4:
  //   identity, 180° rot, 90° CW, 90° CCW, H-flip, V-flip, diagonal,
  //   anti-diagonal
  using T = function<pair<int, int>(int, int)>;
  T transforms[8] = {
      [](int r, int c) { return make_pair(r, c); },
      [R, C](int r, int c) { return make_pair(R - 1 - r, C - 1 - c); },
      [R](int r, int c) { return make_pair(c, R - 1 - r); },
      [C](int r, int c) { return make_pair(C - 1 - c, r); },
      [C](int r, int c) { return make_pair(r, C - 1 - c); },
      [R](int r, int c) { return make_pair(R - 1 - r, c); },
      [](int r, int c) { return make_pair(c, r); },
      [R, C](int r, int c) { return make_pair(C - 1 - c, R - 1 - r); },
  };

  g_numSymm = 0;
  for (int t = 0; t < 8; t++) {
    int perm[49];
    bool valid = true;
    for (int i = 0; i < g_numPos && valid; i++) {
      auto [r, c] = GRID.positions[i];
      auto [nr, nc] = transforms[t](r, c);
      auto it = GRID.coordToIdx.find({nr, nc});
      if (it == GRID.coordToIdx.end())
        valid = false;
      else
        perm[i] = it->second;
    }
    if (valid) {
      for (int i = 0; i < g_numPos; i++)
        g_symmPerm[g_numSymm][i] = perm[i];
      g_numSymm++;
    }
  }
  // (symmetry count intentionally not printed; keep output compact)
}

static void buildLookupTables() {
  g_numPos = GRID.numPositions();
  for (int i = 0; i < g_numPos; i++) {
    const auto &nb = GRID.neighbors[i];
    int count = (int)nb.size();
    g_neighborCount[i] = count;
    g_maxType[i] = GRID.maxType[i];
    for (int j = 0; j < count; j++)
      g_neighborShifts[i][j] = nb[j] * 3;
    for (int j = count; j < 4; j++)
      g_neighborShifts[i][j] = 63; // safe padding: reads always-zero bits
  }
  g_baseState = StateType(0);
  for (int i = 0; i < g_numPos; i++)
    g_baseState = stateSet(g_baseState, i, 1);
}

// ── Path reconstruction helpers ──────────────────────────────────────────────
int scoreGrid(StateType state) {
  int score = 0;
  for (int i = 0; i < g_numPos; i++)
    score += stateGet(state, i);
  return score;
}

// expandPathForEmptyBoard: converts a path from the solver's "all-1 = blue"
// model to a path that starts from a truly empty board.
//
// For each step (pos, newType), before we can upgrade pos, we need:
//   1. pos itself to be explicitly placed (as blue) if not already
//   2. For each required neighbor type 1..newType-1, at least one neighbor
//      must be explicitly placed with that type.
//      - Neighbors with type > 1 were placed by earlier explicit steps.
//      - Neighbors with type == 1 are "implicit blue" and need an explicit
//        "place blue" step inserted here.
//
// Tracks a physical (non-canonical) state so neighbor lookups against
// GRID.neighbors[] remain consistent with position indices.
static vector<pair<int, int>> expandPathForEmptyBoard(const vector<pair<int, int>> &steps) {
  vector<pair<int, int>> out;
  out.reserve(steps.size() * 2);
  uint64_t placed = 0; // bit i = physical position i has been explicitly placed

  // Track both states exactly like printExpandedPath does.
  StateType canonState = g_baseState;
  StateType physState = g_baseState;

  for (auto [pos, newType] : steps) {
    // Find which symmetry maps physState to its canonical form.
    int usedSymm = 0;
    for (int k = 0; k < g_numSymm; k++) {
      if (applyPerm(physState, g_symmPerm[k]) == canonState) {
        usedSymm = k;
        break;
      }
    }
    // perm[i]    = canonical position of physical position i
    // invPerm[j] = physical position of canonical position j
    int invPerm[49];
    for (int i = 0; i < g_numPos; i++)
      invPerm[g_symmPerm[usedSymm][i]] = i;
    int physPos = invPerm[pos];

    // For each required neighbor type req in 1..newType-1, ensure at least one
    // placed physical neighbor has that type. Unplaced neighbors with type 1
    // are "implicit blue"; insert an explicit place-blue step if needed.
    for (int req = 1; req < newType; req++) {
      bool satisfied = false;
      int candidatePhysNb = -1;
      for (int physNb : GRID.neighbors[physPos]) {
        int nbType = stateGet(physState, physNb);
        if (nbType == req) {
          if ((placed >> physNb) & 1) {
            satisfied = true;
            break;
          } else if (req == 1) {
            candidatePhysNb = physNb;
          }
        }
      }
      if (!satisfied && candidatePhysNb != -1) {
        // Express the placement in canonical space for printExpandedPath.
        // Blue placements don't change canonState (already 1), so usedSymm
        // stays stable for this step.
        int canonNb = g_symmPerm[usedSymm][candidatePhysNb];
        out.push_back({canonNb, 1});
        placed |= (1ull << candidatePhysNb);
      }
    }

    out.push_back({pos, newType});
    placed |= (1ull << physPos);
    canonState = canonical(stateSet(canonState, pos, newType));
    physState = stateSet(physState, physPos, newType);
  }

  return out;
}

// result[i] = s[perm[i]] (inverse of applyPerm for display: canon in frame k)
static StateType applyInversePerm(StateType s, const int *perm) {
  StateType result = 0;
  for (int i = 0; i < g_numPos; i++)
    result = stateSet(result, i, stateGet(s, perm[i]));
  return result;
}

static void advanceStateForDisplay(StateType &canonState, StateType &displayState, uint64_t &placedMask, int pos,
                                   int newType) {
  int usedSymm = 0;
  for (int k = 0; k < g_numSymm; k++) {
    if (applyPerm(displayState, g_symmPerm[k]) == canonState) {
      usedSymm = k;
      break;
    }
  }
  int inversePerm[49];
  for (int i = 0; i < g_numPos; i++)
    inversePerm[g_symmPerm[usedSymm][i]] = i;
  int physPos = inversePerm[pos];
  canonState = canonical(stateSet(canonState, pos, newType));
  displayState = stateSet(displayState, physPos, newType);
  placedMask |= (1ull << physPos);
}

// Try all 8 (or g_numSymm) frames; return k whose replay ends at finalCanonState, or -1.
static int chooseFrameMatchingFinal(const vector<pair<int, int>> &expandedSteps, StateType startState,
                                    StateType finalCanonState) {
  for (int k = 0; k < g_numSymm; k++) {
    int invPerm[49];
    for (int i = 0; i < g_numPos; i++)
      invPerm[g_symmPerm[k][i]] = i;
    StateType display = applyInversePerm(startState, g_symmPerm[k]);
    for (const auto &[pos, newType] : expandedSteps) {
      int physPos = invPerm[pos];
      display = stateSet(display, physPos, newType);
    }
    if (display == finalCanonState)
      return k;
  }
  return -1;
}

static void printExpandedPath(const vector<pair<int, int>> &expandedSteps, StateType startState, ostream &out = cout) {
  const bool fromEmpty = (startState == StateType(0));
  StateType canonState = startState;
  StateType displayState = startState;
  uint64_t placedMask = 0;

  // Pre-pass to compute total steps and final canonical state.
  for (auto [pos, newType] : expandedSteps)
    advanceStateForDisplay(canonState, displayState, placedMask, pos, newType);

  StateType finalCanonState = canonState;
  int totalSteps = (int)expandedSteps.size();
  if (!fromEmpty) {
    int missingCount = 0;
    for (int i = 0; i < g_numPos; i++) {
      if (!((placedMask >> i) & 1))
        missingCount++;
    }
    totalSteps += missingCount;
  }

  // Pick the frame whose replay ends at final canonical state (all 8 are cheap).
  int fixedSymm = chooseFrameMatchingFinal(expandedSteps, startState, finalCanonState);
  int inversePerm[49];
  if (fixedSymm >= 0) {
    for (int i = 0; i < g_numPos; i++)
      inversePerm[g_symmPerm[fixedSymm][i]] = i;
  }

  // Reset for actual print pass.
  canonState = startState;
  displayState = (fixedSymm >= 0) ? applyInversePerm(startState, g_symmPerm[fixedSymm]) : startState;
  placedMask = 0;

  out << "\n--- Build path (" << totalSteps << " steps) ---\n\n";
  out << "Start (" << (fromEmpty ? "empty" : "all blue") << " board):\n";
  if (fixedSymm >= 0)
    printGrid(applyInversePerm(startState, g_symmPerm[fixedSymm]), uint64_t(0), out);
  else
    printGrid(startState, uint64_t(0), out);
  int validSteps = 0, invalidSteps = 0;
  int stepNum = 1;

  for (auto [pos, newType] : expandedSteps) {
    int usedSymm;
    int physPos;
    if (fixedSymm >= 0) {
      usedSymm = fixedSymm;
      physPos = inversePerm[pos];
    } else {
      int k = 0;
      for (; k < g_numSymm; k++) {
        if (applyPerm(displayState, g_symmPerm[k]) == canonState)
          break;
      }
      usedSymm = k;
      for (int i = 0; i < g_numPos; i++)
        inversePerm[g_symmPerm[usedSymm][i]] = i;
      physPos = inversePerm[pos];
    }
    bool wasPlaced = (placedMask >> physPos) & 1;
    int oldType = wasPlaced ? stateGet(displayState, physPos) : -1;

    // Validity check in canonical space
    uint32_t present = neighborPresent(canonState, pos);
    bool valid = canPlace(present, newType);
    if (valid)
      validSteps++;
    else
      invalidSteps++;

    advanceStateForDisplay(canonState, displayState, placedMask, pos, newType);

    auto [r, c] = GRID.positions[physPos];
    out << "\nStep " << stepNum << ": pos " << physPos << " (" << r << "," << c << ")  ";
    if (oldType == -1) {
      out << "empty -> " << COLORS[newType] << BLACK_TEXT << COLOR_NAMES[newType] << RESET << " (" << newType << ")";
    } else {
      out << COLORS[oldType] << BLACK_TEXT << COLOR_NAMES[oldType] << RESET << " (" << oldType << ") -> "
          << COLORS[newType] << BLACK_TEXT << COLOR_NAMES[newType] << RESET << " (" << newType << ")";
    }
    out << (valid ? "  \033[32m\u2713\033[0m" : "  \033[31m\u2717 INVALID\033[0m") << "\n";
    printGrid(displayState, placedMask, out);
    stepNum++;
  }

  // When path started from all-blue: any position never explicitly placed is an
  // implicit blue; show a final step for each.
  if (!fromEmpty) {
    for (int i = 0; i < g_numPos; i++) {
      if (!((placedMask >> i) & 1)) {
        placedMask |= (1ull << i);
        auto [r, c] = GRID.positions[i];
        out << "\nStep " << stepNum << ": pos " << i << " (" << r << "," << c << ")  ";
        out << "empty -> " << COLORS[1] << BLACK_TEXT << COLOR_NAMES[1] << RESET << " (1)";
        out << "  \033[32m\u2713\033[0m"
            << "\n";
        printGrid(displayState, placedMask, out);
        stepNum++;
        validSteps++;
      }
    }
  }

  int reachedScore = scoreGrid(displayState);
  out << "\nBuild path: " << totalSteps << " steps, " << validSteps << " valid, " << invalidSteps << " invalid";
  if (invalidSteps == 0) {
    out << " \033[32m\u2014 ALL VALID\033[0m";
  } else {
    out << " \033[31m\u2014 HAS INVALID STEPS\033[0m";
  }
  out << "\nReached score: " << reachedScore << "\n";
}

template <typename S>
inline void printBuildPath(StateType bestState,
                           const ankerl::unordered_dense::segmented_set<S, StorageHash<S>, StorageEqual<S>> &visited,
                           ostream &out = cout) {
  vector<pair<int, int>> steps; // (pos, newType) — pos in state *before* the move
  StateType cur = bestState;
  while (true) {
    auto tr = storageTrans<S>(*visited.find(toStorageKey<S>(cur)));
    if (tr.getPos() == -1)
      break;
    int childPos = tr.getPos();
    int newType = stateGet(cur, childPos);
    StateType parent = canonical(stateSet(cur, childPos, tr.getType()));
    // Replay expects (pos, newType) where pos is in the state before the move.
    // childPos is in cur (child); find pos in parent s.t. stateSet(parent, pos,
    // newType) -> cur.
    int posInParent = -1;
    for (int p = 0; p < g_numPos; p++) {
      if (stateGet(parent, p) == newType)
        continue;
      if (canonical(stateSet(parent, p, newType)) == cur) {
        posInParent = p;
        break;
      }
    }
    if (posInParent < 0)
      posInParent = childPos; // fallback if no symmetry (g_numSymm == 1)
    steps.push_back({posInParent, newType});
    cur = parent;
  }
  reverse(steps.begin(), steps.end());
  if (PATH_FROM_EMPTY)
    printExpandedPath(steps, StateType(0), out);
  else
    printExpandedPath(expandPathForEmptyBoard(steps), g_baseState, out);
}

// Return path (steps) from start to target using the BFS visited set. Empty if target not in visited.
template <typename S>
static vector<pair<int, int>>
getPathFromVisited(StateType target,
                   const ankerl::unordered_dense::segmented_set<S, StorageHash<S>, StorageEqual<S>> &visited) {
  vector<pair<int, int>> steps;
  StateType cur = target;
  while (true) {
    auto it = visited.find(toStorageKey<S>(cur));
    if (it == visited.end())
      return {};
    auto tr = storageTrans<S>(*it);
    if (tr.getPos() == -1)
      break;
    int childPos = tr.getPos();
    int newType = stateGet(cur, childPos);
    StateType parent = canonical(stateSet(cur, childPos, tr.getType()));
    int posInParent = -1;
    for (int p = 0; p < g_numPos; p++) {
      if (stateGet(parent, p) == newType)
        continue;
      if (canonical(stateSet(parent, p, newType)) == cur) {
        posInParent = p;
        break;
      }
    }
    if (posInParent < 0)
      posInParent = childPos;
    steps.push_back({posInParent, newType});
    cur = parent;
  }
  reverse(steps.begin(), steps.end());
  return steps;
}

void printBuildPathSteps(const vector<pair<int, int>> &buildPath, ostream &out) {
  if (PATH_FROM_EMPTY)
    printExpandedPath(buildPath, StateType(0), out);
  else
    printExpandedPath(expandPathForEmptyBoard(buildPath), g_baseState, out);
}

// BFS from base state (all blue) stopping as soon as target is reached.
// Returns a set of encoded states (transition packed into top byte) for path
// reconstruction.
inline ankerl::unordered_dense::segmented_set<StateType, StateHash, StateEqual> bfsShortestPath(StateType target) {
  ankerl::unordered_dense::segmented_set<StateType, StateHash, StateEqual> path;
  struct Entry {
    StateType state;
    int score;
  };
  vector<Entry> q;
  size_t head = 0;

  StateType bfsStart = PATH_FROM_EMPTY ? StateType(0) : g_baseState;
  path.emplace(encodeTransition(bfsStart, {})); // sentinel: packed=255
  q.push_back({bfsStart, scoreGrid(bfsStart)});

  while (head < q.size()) {
    auto [cur, curScore] = q[head++];
    if (cur == target)
      break;
    for (int pos = 0; pos < g_numPos; pos++) {
      int curType = stateGet(cur, pos);
      uint32_t present = neighborPresent(cur, pos);
      for (int newType = 0; newType <= g_maxType[pos]; newType++) {
        if (!FORWARD_TRANSITIONS[curType][newType])
          continue;
        if (!canPlace(present, newType))
          continue;
        auto [next, childPos] = canonicalAndChildPos(stateSet(cur, pos, newType), pos);
        if (path.emplace(encodeTransition(next, Transition{(int8_t)childPos, (int8_t)curType})).second)
          q.push_back({next, curScore + (newType - curType)});
      }
    }
  }
  return path;
}

// ── Best-states file (read/write) ─────────────────────────────────────────────
// Line format: <state_digits> <path_length> <state_0> <state_1> ... <state_n>
// (path_length = number of moves; state_0..state_n = full sequence from start to goal.)
// Legacy: line may be state_digits only (no path); we then run reverse BFS to fill path.

struct BestStateEntry {
  string stateStr;              // g_numPos digits
  int pathLen = -1;             // -1 = missing
  vector<string> stateSequence; // pathLen+1 state strings when pathLen >= 0
};

// Serialize one canonical state to a string of g_numPos digits '0'-'5'.
static string stateToLine(StateType s) {
  string line;
  line.reserve((size_t)g_numPos);
  for (int pos = 0; pos < g_numPos; pos++)
    line += (char)('0' + stateGet(s, pos));
  return line;
}

// Parse first g_numPos digits of s as canonical state.
static StateType lineToState(const string &s) {
  StateType st = 0;
  for (int i = 0; i < g_numPos && i < (int)s.size(); i++) {
    int t = s[i] - '0';
    if (t < 0 || t > 5)
      return 0;
    st = stateSet(st, i, t);
  }
  return canonical(st);
}

// From startState apply each step; return sequence of states (startState, ..., endState).
static vector<StateType> stateSequenceFromSteps(StateType startState, const vector<pair<int, int>> &steps) {
  vector<StateType> seq;
  seq.reserve(steps.size() + 1);
  seq.push_back(startState);
  StateType cur = startState;
  for (const auto &[pos, newType] : steps) {
    cur = canonical(stateSet(cur, pos, newType));
    seq.push_back(cur);
  }
  return seq;
}

// Read best_states file. Returns (score, entries). Invalid/missing -> (-1, {}).
static pair<int, vector<BestStateEntry>> readBestStatesFile(const string &filename) {
  ifstream f(filename);
  if (!f)
    return {-1, {}};
  string header;
  if (!getline(f, header) || header.find("score ") != 0)
    return {-1, {}};
  int score = 0;
  if (sscanf(header.c_str(), "score %d", &score) != 1)
    return {-1, {}};
  string countLine;
  if (!getline(f, countLine) || countLine.find("count ") != 0)
    return {-1, {}};
  int count = 0;
  if (sscanf(countLine.c_str(), "count %d", &count) != 1 || count < 0)
    return {-1, {}};
  vector<BestStateEntry> entries;
  entries.reserve((size_t)min(count, 10000));
  for (int i = 0; i < count; i++) {
    string line;
    if (!getline(f, line))
      return {-1, {}};
    BestStateEntry e;
    size_t pos = 0;
    while (pos < line.size() && line[pos] == ' ')
      pos++;
    if (pos >= line.size() || (int)line.size() - (int)pos < g_numPos)
      return {-1, {}};
    e.stateStr = line.substr(pos, (size_t)g_numPos);
    pos += (size_t)g_numPos;
    for (int j = 0; j < g_numPos && j < (int)e.stateStr.size(); j++)
      if (e.stateStr[j] < '0' || e.stateStr[j] > '5')
        return {-1, {}};
    while (pos < line.size() && line[pos] == ' ')
      pos++;
    if (pos >= line.size()) {
      e.pathLen = -1;
      entries.push_back(std::move(e));
      continue;
    }
    // Parse path_length and state sequence
    size_t tokStart = pos;
    while (pos < line.size() && line[pos] != ' ')
      pos++;
    string pathLenStr = line.substr(tokStart, pos - tokStart);
    int pathLen = 0;
    if (sscanf(pathLenStr.c_str(), "%d", &pathLen) != 1 || pathLen < 0)
      return {-1, {}};
    e.pathLen = pathLen;
    while (pos < line.size() && line[pos] == ' ')
      pos++;
    for (int k = 0; k <= pathLen && pos < line.size(); k++) {
      if (pos + (size_t)g_numPos > line.size())
        return {-1, {}};
      string stateTok = line.substr(pos, (size_t)g_numPos);
      for (int j = 0; j < g_numPos; j++)
        if (stateTok[j] < '0' || stateTok[j] > '5')
          return {-1, {}};
      e.stateSequence.push_back(stateTok);
      pos += (size_t)g_numPos;
      while (pos < line.size() && line[pos] == ' ')
        pos++;
    }
    if ((int)e.stateSequence.size() != pathLen + 1)
      return {-1, {}};
    entries.push_back(std::move(e));
  }
  return {score, entries};
}

// ── Sliding Window Optimization ─────────────────────────────────────────────-
// Two window placement policies are supported:
// - Center-based (original): choose a window *center* at an existing grid cell, require the kernel's
//   bounding box to fit inside the grid rectangle, and build windows by centering the pattern.
// - Top-left-based (new, rectangular kernels): choose a window *top-left* anywhere in the grid rectangle,
//   skip placements based on per-strip overlap against the actual level pattern, and build windows
//   by placing the pattern starting at (r0,c0).

struct KernelSize {
  int rows = 0;
  int cols = 0;
};

static KernelSize getKernelSize() {
  const auto &pat = WINDOW_KERNEL_PATTERN;
  KernelSize k;
  k.rows = (int)pat.size();
  for (const auto &row : pat)
    k.cols = max(k.cols, (int)row.size());
  return k;
}

// Rotate the active kernel pattern 90 degrees clockwise. Computed once and cached.
static const vector<string> &getWindowKernelPatternRotated90() {
  static vector<string> rotated;
  if (!rotated.empty())
    return rotated;

  const auto &pat = WINDOW_KERNEL_PATTERN;
  int R = (int)pat.size();
  int C = 0;
  for (const auto &row : pat)
    C = max(C, (int)row.size());

  rotated.assign(C, string(R, ' '));
  for (int r = 0; r < R; r++) {
    for (int c = 0; c < (int)pat[r].size(); c++) {
      char ch = pat[r][c];
      int nr = c;
      int nc = R - 1 - r;
      rotated[nr][nc] = ch;
    }
  }
  return rotated;
}

// Shared helper: redraw the in-place "New best" block (score + grid).
static void redrawNewBestBlock(int bestScore, StateType bestState) {
  if (g_activeStatusLine)
    g_activeStatusLine->clear();
  vector<string> lines;
  lines.push_back("New best: " + to_string(bestScore));
  vector<string> gridLines = gridToLines(bestState);
  for (auto &l : gridLines)
    lines.push_back(l);
  if (s_newBestBlockLines > 0)
    cout << "\033[" << s_newBestBlockLines << "A";
  for (size_t i = 0; i < lines.size(); i++)
    cout << "\r\033[K" << lines[i] << "\n";
  cout << flush;
  s_newBestBlockLines = (int)lines.size();
  if (g_activeStatusLine)
    g_activeStatusLine->redraw();
}

// Shared BFS core over a fixed set of window positions.
// When bestCanonicalStates is non-null, collect up to maxBestStates canonical states that achieve the best score.
static void
exploreWindowBfs(StateType startState, const vector<int> &windowPositions, vector<StateType> &newStates, int &bestScore,
                 StateType &bestState, long long &explored,
                 ankerl::unordered_dense::set<StateType, StateHash, StateEqual> *bestCanonicalStates = nullptr,
                 int maxBestStates = SLIDING_WINDOW_MAX_BEST_STATES) {
  ankerl::unordered_dense::segmented_set<StateType, StateHash, StateEqual> localVisited;
  struct Entry {
    StateType state;
    int score;
    int depth;
  };
  vector<Entry> q;

  localVisited.insert(startState);
  q.push_back({startState, scoreGrid(startState), 0});
  size_t head = 0;

  while (head < q.size()) {
    auto [cur, curScore, depth] = q[head++];
    explored++;

    if (depth >= MAX_DEPTH)
      continue;

    if (curScore > bestScore) {
      bestScore = curScore;
      bestState = cur;
      if (bestCanonicalStates) {
        bestCanonicalStates->clear();
        bestCanonicalStates->insert(canonical(cur));
      }
      redrawNewBestBlock(bestScore, bestState);
    } else if (curScore == bestScore && bestCanonicalStates && (int)bestCanonicalStates->size() < maxBestStates) {
      bestCanonicalStates->insert(canonical(cur));
    }

    for (int pos : windowPositions) {
      if (pos >= g_numPos)
        continue;
      int curType = stateGet(cur, pos);

      uint32_t present = neighborPresent(cur, pos);
      for (int newType = 0; newType <= g_maxType[pos]; newType++) {
        if (!FORWARD_TRANSITIONS_SLIDING_WINDOW[curType][newType])
          continue;
        if (!canPlace(present, newType))
          continue;

        StateType next = stateSet(cur, pos, newType);
        if (localVisited.insert(next).second) {
          int nextScore = curScore + (newType - curType);
          q.push_back({next, nextScore, depth + 1});
          newStates.push_back(next);
        }
      }
    }
  }
}

// ── Original center-based path/window (used for non-rectangular kernels) ──

// Build zigzag path of *center indices* starting from a given corner (0..3).
static vector<int> buildZigzagPathCenter(int corner = 0) {
  const int R = GRID.rows, C = GRID.cols;
  vector<vector<pair<int, int>>> diagonals;

  for (int r = 0; r < R; r++) {
    for (int c = 0; c < C; c++) {
      if (!GRID.coordToIdx.count({r, c}))
        continue; // center must be a real grid cell

      int diag;
      switch (corner) {
      case 0:
        diag = r + c;
        break; // top-left
      case 1:
        diag = r + (C - 1 - c);
        break; // top-right
      case 2:
        diag = (R - 1 - r) + c;
        break; // bottom-left
      case 3:
        diag = (R - 1 - r) + (C - 1 - c);
        break; // bottom-right
      default:
        diag = r + c;
        break;
      }
      if (diag >= (int)diagonals.size())
        diagonals.resize(diag + 1);
      diagonals[diag].push_back({r, c});
    }
  }

  vector<int> path;
  for (int d = 0; d < (int)diagonals.size(); d++) {
    auto &diag = diagonals[d];
    if (d % 2 == 1)
      reverse(diag.begin(), diag.end());
    for (auto [r, c] : diag) {
      path.push_back(GRID.coordToIdx.at({r, c}));
    }
  }
  return path;
}

// Build snake path of *center indices*: row 0 L→R, row 1 R→L, ...
static vector<int> buildSnakePathCenter() {
  const int R = GRID.rows, C = GRID.cols;
  vector<vector<int>> rowCols(R);
  for (int r = 0; r < R; r++) {
    for (int c = 0; c < C; c++) {
      if (!GRID.coordToIdx.count({r, c}))
        continue;
      rowCols[r].push_back(c);
    }
    if (r % 2 == 1)
      reverse(rowCols[r].begin(), rowCols[r].end());
  }
  vector<int> path;
  for (int r = 0; r < R; r++) {
    for (int c : rowCols[r])
      path.push_back(GRID.coordToIdx.at({r, c}));
  }
  return path;
}

// Center-based window exploration: centerIdx must be a valid grid position index.
static void
exploreWindowCenter(StateType startState, int centerIdx, vector<StateType> &newStates, int &bestScore,
                    StateType &bestState, long long &explored,
                    ankerl::unordered_dense::set<StateType, StateHash, StateEqual> *bestCanonicalStates = nullptr,
                    int maxBestStates = SLIDING_WINDOW_MAX_BEST_STATES) {
  auto [cr, cc] = GRID.positions[centerIdx];
  const auto &pat = WINDOW_KERNEL_PATTERN;
  KernelSize k = getKernelSize();
  int centerR = k.rows / 2;
  int centerC = k.cols / 2;

  vector<int> windowPositions;
  for (int pr = 0; pr < (int)pat.size(); pr++) {
    for (int pc = 0; pc < (int)pat[pr].size(); pc++) {
      if (pat[pr][pc] != '.')
        continue;
      int dr = pr - centerR;
      int dc = pc - centerC;
      auto it = GRID.coordToIdx.find({cr + dr, cc + dc});
      if (it != GRID.coordToIdx.end())
        windowPositions.push_back(it->second);
    }
  }

  printGridWithWindow(windowPositions);
  cout << "\n" << flush;
  exploreWindowBfs(startState, windowPositions, newStates, bestScore, bestState, explored, bestCanonicalStates,
                   maxBestStates);
}

// Bounding-box fit: require the kernel's rectangular extent from (topLeftRow, topLeftCol)
// to (topLeftRow + rows - 1, topLeftCol + cols - 1) to lie inside the grid rectangle.
static bool shouldSkipTopLeftBoundingBox(int topLeftRow, int topLeftCol, const KernelSize &k) {
  const int R = GRID.rows, C = GRID.cols;
  int rLo = topLeftRow;
  int rHi = topLeftRow + k.rows - 1;
  int cLo = topLeftCol;
  int cHi = topLeftCol + k.cols - 1;
  if (rLo < 0 || cLo < 0 || rHi >= R || cHi >= C)
    return true;
  return false;
}

// Pattern-overlap fit: for a given top-left, a strip "overlaps the grid" iff at least one '.'
// cell in that strip lands on a real grid cell (per createGridFromPattern()).
// Skip this top-left if any dotted strip has zero overlap.
static bool shouldSkipTopLeftPatternOverlapWithPattern(int topLeftRow, int topLeftCol, const KernelSize &k,
                                                       const vector<string> &pat) {
  for (int kr = 0; kr < k.rows; kr++) {
    bool rowHasDot = false;
    bool rowOverlaps = false;
    if (kr < (int)pat.size()) {
      for (int pc = 0; pc < (int)pat[kr].size(); pc++) {
        if (pat[kr][pc] != '.')
          continue;
        rowHasDot = true;
        if (GRID.coordToIdx.count({topLeftRow + kr, topLeftCol + pc})) {
          rowOverlaps = true;
          break;
        }
      }
    }
    if (rowHasDot && !rowOverlaps)
      return true;
  }
  for (int kc = 0; kc < k.cols; kc++) {
    bool colHasDot = false;
    bool colOverlaps = false;
    for (int kr = 0; kr < (int)pat.size(); kr++) {
      if (kc >= (int)pat[kr].size())
        continue;
      if (pat[kr][kc] != '.')
        continue;
      colHasDot = true;
      if (GRID.coordToIdx.count({topLeftRow + kr, topLeftCol + kc})) {
        colOverlaps = true;
        break;
      }
    }
    if (colHasDot && !colOverlaps)
      return true;
  }
  return false;
}

static bool shouldSkipTopLeftPatternOverlap(int topLeftRow, int topLeftCol, const KernelSize &k) {
  const auto &pat = WINDOW_KERNEL_PATTERN;
  return shouldSkipTopLeftPatternOverlapWithPattern(topLeftRow, topLeftCol, k, pat);
}

// Dispatch between fit modes based on kernel index (see getWindowKernelFitMode()).
static bool shouldSkipTopLeft(int topLeftRow, int topLeftCol, const KernelSize &k) {
  KernelFitMode mode = getWindowKernelFitMode();
  if (mode == KernelFitMode::BoundingBox)
    return shouldSkipTopLeftBoundingBox(topLeftRow, topLeftCol, k);
  return shouldSkipTopLeftPatternOverlap(topLeftRow, topLeftCol, k);
}

// Collect every (r, c) in the grid rectangle that is not skipped for the active kernel.
static vector<pair<int, int>> collectValidTopLefts(const KernelSize &k) {
  const int R = GRID.rows, C = GRID.cols;
  vector<pair<int, int>> out;
  for (int r = 0; r < R; r++) {
    for (int c = 0; c < C; c++) {
      if (shouldSkipTopLeft(r, c, k))
        continue;
      out.push_back({r, c});
    }
  }
  return out;
}

// Collect every (r, c) in the grid rectangle that is not skipped for a specific pattern (pattern-overlap mode).
static vector<pair<int, int>> collectValidTopLeftsForPattern(const KernelSize &k, const vector<string> &pat) {
  const int R = GRID.rows, C = GRID.cols;
  vector<pair<int, int>> out;
  for (int r = 0; r < R; r++) {
    for (int c = 0; c < C; c++) {
      if (shouldSkipTopLeftPatternOverlapWithPattern(r, c, k, pat))
        continue;
      out.push_back({r, c});
    }
  }
  return out;
}

// Order by diagonal (corner 0..3), zigzag within each diagonal; return (row, col) top-lefts.
static vector<pair<int, int>> buildZigzagPathTopLeft(int corner = 0) {
  const int R = GRID.rows, C = GRID.cols;
  KernelSize k = getKernelSize();
  vector<pair<int, int>> cells = collectValidTopLefts(k);

  auto diagIndex = [&](int r, int c) {
    switch (corner) {
    case 0:
      return r + c;
    case 1:
      return r + (C - 1 - c);
    case 2:
      return (R - 1 - r) + c;
    case 3:
      return (R - 1 - r) + (C - 1 - c);
    default:
      return r + c;
    }
  };

  sort(cells.begin(), cells.end(), [&](const pair<int, int> &a, const pair<int, int> &b) {
    int da = diagIndex(a.first, a.second), db = diagIndex(b.first, b.second);
    if (da != db)
      return da < db;
    // Zigzag within diagonal: alternate direction by diagonal index
    if (da % 2 == 0)
      return a.second < b.second;
    return a.second > b.second;
  });
  return cells;
}

static vector<pair<int, int>> buildZigzagPathTopLeftForPattern(const vector<string> &pat, int corner = 0) {
  const int R = GRID.rows, C = GRID.cols;
  KernelSize k;
  k.rows = (int)pat.size();
  k.cols = 0;
  for (const auto &row : pat)
    k.cols = max(k.cols, (int)row.size());
  vector<pair<int, int>> cells = collectValidTopLeftsForPattern(k, pat);

  auto diagIndex = [&](int r, int c) {
    switch (corner) {
    case 0:
      return r + c;
    case 1:
      return r + (C - 1 - c);
    case 2:
      return (R - 1 - r) + c;
    case 3:
      return (R - 1 - r) + (C - 1 - c);
    default:
      return r + c;
    }
  };

  sort(cells.begin(), cells.end(), [&](const pair<int, int> &a, const pair<int, int> &b) {
    int da = diagIndex(a.first, a.second), db = diagIndex(b.first, b.second);
    if (da != db)
      return da < db;
    if (da % 2 == 0)
      return a.second < b.second;
    return a.second > b.second;
  });
  return cells;
}

// Order by row; within each row, left-to-right or right-to-left (snake); return (row, col) top-lefts.
static vector<pair<int, int>> buildSnakePathTopLeft() {
  KernelSize k = getKernelSize();
  vector<pair<int, int>> cells = collectValidTopLefts(k);

  sort(cells.begin(), cells.end(), [](const pair<int, int> &a, const pair<int, int> &b) {
    if (a.first != b.first)
      return a.first < b.first;
    return a.second < b.second;
  });

  for (int r = 0; r < GRID.rows; r++) {
    auto rowStart = lower_bound(cells.begin(), cells.end(), make_pair(r, 0));
    auto rowEnd = lower_bound(cells.begin(), cells.end(), make_pair(r + 1, 0));
    if (r % 2 == 1 && rowStart != rowEnd)
      reverse(rowStart, rowEnd);
  }
  return cells;
}

static vector<pair<int, int>> buildSnakePathTopLeftForPattern(const vector<string> &pat) {
  KernelSize k;
  k.rows = (int)pat.size();
  k.cols = 0;
  for (const auto &row : pat)
    k.cols = max(k.cols, (int)row.size());
  vector<pair<int, int>> cells = collectValidTopLeftsForPattern(k, pat);

  sort(cells.begin(), cells.end(), [](const pair<int, int> &a, const pair<int, int> &b) {
    if (a.first != b.first)
      return a.first < b.first;
    return a.second < b.second;
  });

  for (int r = 0; r < GRID.rows; r++) {
    auto rowStart = lower_bound(cells.begin(), cells.end(), make_pair(r, 0));
    auto rowEnd = lower_bound(cells.begin(), cells.end(), make_pair(r + 1, 0));
    if (r % 2 == 1 && rowStart != rowEnd)
      reverse(rowStart, rowEnd);
  }
  return cells;
}

// BFS within a window. (topLeftRow, topLeftCol) is the kernel's top-left; it may be a hole.
static void
exploreWindowTopLeft(StateType startState, int topLeftRow, int topLeftCol, vector<StateType> &newStates, int &bestScore,
                     StateType &bestState, long long &explored,
                     ankerl::unordered_dense::set<StateType, StateHash, StateEqual> *bestCanonicalStates = nullptr,
                     int maxBestStates = SLIDING_WINDOW_MAX_BEST_STATES) {
  int r0 = topLeftRow, c0 = topLeftCol;
  KernelSize k = getKernelSize();
  const auto &pat = WINDOW_KERNEL_PATTERN;

  vector<int> windowPositions;
  for (int pr = 0; pr < k.rows; pr++) {
    for (int pc = 0; pc < (int)pat[pr].size(); pc++) {
      if (pat[pr][pc] != '.')
        continue;
      auto it = GRID.coordToIdx.find({r0 + pr, c0 + pc});
      if (it != GRID.coordToIdx.end())
        windowPositions.push_back(it->second);
    }
  }

  printGridWithWindow(windowPositions);
  cout << "\n" << flush;

  exploreWindowBfs(startState, windowPositions, newStates, bestScore, bestState, explored, bestCanonicalStates,
                   maxBestStates);
}

// Top-left exploration with the kernel rotated 90 degrees clockwise.
static void exploreWindowTopLeftRotated(
    StateType startState, int topLeftRow, int topLeftCol, vector<StateType> &newStates, int &bestScore,
    StateType &bestState, long long &explored,
    ankerl::unordered_dense::set<StateType, StateHash, StateEqual> *bestCanonicalStates = nullptr,
    int maxBestStates = SLIDING_WINDOW_MAX_BEST_STATES) {
  int r0 = topLeftRow, c0 = topLeftCol;
  const auto &pat = getWindowKernelPatternRotated90();
  KernelSize k;
  k.rows = (int)pat.size();
  k.cols = 0;
  for (const auto &row : pat)
    k.cols = max(k.cols, (int)row.size());

  vector<int> windowPositions;
  for (int pr = 0; pr < k.rows; pr++) {
    for (int pc = 0; pc < (int)pat[pr].size(); pc++) {
      if (pat[pr][pc] != '.')
        continue;
      auto it = GRID.coordToIdx.find({r0 + pr, c0 + pc});
      if (it != GRID.coordToIdx.end())
        windowPositions.push_back(it->second);
    }
  }

  printGridWithWindow(windowPositions);
  cout << "\n" << flush;

  exploreWindowBfs(startState, windowPositions, newStates, bestScore, bestState, explored, bestCanonicalStates,
                   maxBestStates);
}

static SearchResult runSlidingWindow() {
  s_newBestBlockLines = 0;
  int baseScore = scoreGrid(g_baseState);
  SearchResult result{g_baseState, baseScore, 0, 0};
  result.bestState = g_baseState;
  result.bestScore = baseScore;
  result.explored = 0;

  // (kernel pattern intentionally not printed; keep output compact)

  const char *cornerNames[] = {"top-left", "top-right", "bottom-left", "bottom-right"};

  struct WindowStats {
    int posIdx;
    int row, col;
    int kernelSize;
    long long statesExplored;
    long long timeMs;
    bool improved;
    int scoreAfter;
  };

  StateType currentState = g_baseState;
  int currentScore = baseScore;
  int roundNum = 0;
  ankerl::unordered_dense::set<StateType, StateHash, StateEqual> bestCanonicalStates;

  string prefix = "";
  if (std::filesystem::current_path().filename() == "build")
    prefix = "../";
  string bestStatesFilename = prefix + "best_states/level" + to_string(getCurrentLevel()) + "_type" +
                              to_string(getCurrentUnlockedType()) + ".txt";
  std::filesystem::create_directories(std::filesystem::path(bestStatesFilename).parent_path());
  int fileBestScore = -1;
  vector<BestStateEntry> fileEntries;
  {
    auto [score, entries] = readBestStatesFile(bestStatesFilename);
    fileBestScore = score;
    fileEntries = std::move(entries);
  }

  while (true) {
    roundNum++;
    int scoreAtRoundStart = currentScore;

    // (round header intentionally not printed; keep output compact)

    int corner = 1; // top-right (used only when WINDOW_PATH_TYPE == 0)
    KernelFitMode fitMode = getWindowKernelFitMode();

    int bestScore = currentScore;
    StateType bestState = currentState;
    long long explored = 0;

    int totalWindows = 0;
    // (pass header intentionally not printed; keep output compact)

    auto passStart = chrono::steady_clock::now();

    vector<WindowStats> stats;

    if (fitMode == KernelFitMode::BoundingBox) {
      // Original behavior: choose window centers at real grid cells.
      vector<int> windowPath = (WINDOW_PATH_TYPE == 1) ? buildSnakePathCenter() : buildZigzagPathCenter(corner);
      totalWindows = (int)windowPath.size();

      for (int windowIdx = 0; windowIdx < totalWindows; windowIdx++) {
        s_newBestBlockLines = 0; // new window: don't cursor-up
        vector<StateType> exploredStates;
        long long exploredBefore = explored;

        int posIdx = windowPath[windowIdx];
        auto [cr, cc] = GRID.positions[posIdx];

        auto t0 = chrono::steady_clock::now();
        exploreWindowCenter(currentState, posIdx, exploredStates, bestScore, bestState, explored, &bestCanonicalStates,
                            SLIDING_WINDOW_MAX_BEST_STATES);
        auto t1 = chrono::steady_clock::now();
        auto ms = chrono::duration_cast<chrono::milliseconds>(t1 - t0).count();

        bool improved = (bestState != currentState);
        if (improved)
          currentState = bestState;

        stats.push_back(
            {posIdx, cr, cc, (int)exploredStates.size(), explored - exploredBefore, ms, improved, bestScore});
      }
    } else {
      // Rectangular-kernel behavior: choose window top-lefts in the grid rectangle,
      // and run two passes: original kernel, then the kernel rotated 90 degrees.
      const auto &pat0 = WINDOW_KERNEL_PATTERN;
      const auto &pat1 = getWindowKernelPatternRotated90();

      vector<pair<int, int>> windowPath0 = (WINDOW_PATH_TYPE == 1) ? buildSnakePathTopLeftForPattern(pat0)
                                                                   : buildZigzagPathTopLeftForPattern(pat0, corner);
      vector<pair<int, int>> windowPath1 = (WINDOW_PATH_TYPE == 1) ? buildSnakePathTopLeftForPattern(pat1)
                                                                   : buildZigzagPathTopLeftForPattern(pat1, corner);

      totalWindows = (int)windowPath0.size() + (int)windowPath1.size();

      // Pass 0: original orientation
      for (const auto &[r0, c0] : windowPath0) {
        s_newBestBlockLines = 0; // new window: don't cursor-up
        vector<StateType> exploredStates;
        long long exploredBefore = explored;

        int posIdx = -1;
        auto itIdx = GRID.coordToIdx.find({r0, c0});
        if (itIdx != GRID.coordToIdx.end())
          posIdx = itIdx->second;

        auto t0 = chrono::steady_clock::now();
        exploreWindowTopLeft(currentState, r0, c0, exploredStates, bestScore, bestState, explored, &bestCanonicalStates,
                             SLIDING_WINDOW_MAX_BEST_STATES);
        auto t1 = chrono::steady_clock::now();
        auto ms = chrono::duration_cast<chrono::milliseconds>(t1 - t0).count();

        bool improved = (bestState != currentState);
        if (improved)
          currentState = bestState;

        stats.push_back(
            {posIdx, r0, c0, (int)exploredStates.size(), explored - exploredBefore, ms, improved, bestScore});
      }

      // Pass 1: rotated 90 degrees
      for (const auto &[r0, c0] : windowPath1) {
        s_newBestBlockLines = 0; // new window: don't cursor-up
        vector<StateType> exploredStates;
        long long exploredBefore = explored;

        int posIdx = -1;
        auto itIdx = GRID.coordToIdx.find({r0, c0});
        if (itIdx != GRID.coordToIdx.end())
          posIdx = itIdx->second;

        auto t0 = chrono::steady_clock::now();
        exploreWindowTopLeftRotated(currentState, r0, c0, exploredStates, bestScore, bestState, explored,
                                    &bestCanonicalStates, SLIDING_WINDOW_MAX_BEST_STATES);
        auto t1 = chrono::steady_clock::now();
        auto ms = chrono::duration_cast<chrono::milliseconds>(t1 - t0).count();

        bool improved = (bestState != currentState);
        if (improved)
          currentState = bestState;

        stats.push_back(
            {posIdx, r0, c0, (int)exploredStates.size(), explored - exploredBefore, ms, improved, bestScore});
      }
    }

    currentScore = bestScore;

    auto passEnd = chrono::steady_clock::now();
    auto passMs = chrono::duration_cast<chrono::milliseconds>(passEnd - passStart).count();

    // Print summary table (built from widths; no hardcoded border strings)
    cout << "\n";
    const vector<int> widths = {5, 9, 14, 16, 12, 6, 7};
    const vector<TableAlign> headerAlign = {TableAlign::Center, TableAlign::Center, TableAlign::Center,
                                            TableAlign::Center, TableAlign::Center, TableAlign::Center,
                                            TableAlign::Center};
    const vector<TableAlign> rowAlign = {TableAlign::Right, TableAlign::Center, TableAlign::Right, TableAlign::Right,
                                         TableAlign::Right, TableAlign::Center, TableAlign::Right};

    printTableBorderTop(widths);
    printTableRow({"Pos", "(r,c)", "Kernel", "States", "Time", "Imp", "Score"}, widths, headerAlign);
    printTableBorderMid(widths);
    for (const auto &s : stats) {
      printTableRow({to_string(s.posIdx), "(" + to_string(s.row) + ", " + to_string(s.col) + ")",
                     formatNum((size_t)s.kernelSize), formatNum((long long)s.statesExplored), formatTimeMs(s.timeMs),
                     (s.improved ? "Y" : "N"), to_string(s.scoreAfter)},
                    widths, rowAlign);
    }
    printTableBorderBottom(widths);
    cout << "Pass " << (char)('A' + corner) << " (" << cornerNames[corner] << "): score " << currentScore << " ("
         << formatTimeMs((long long)passMs) << ", " << formatNum(explored) << " states)\n\n"
         << flush;

    result.explored += explored;

    cout << "Round " << roundNum << " complete: score " << scoreAtRoundStart << " -> " << currentScore << "\n" << flush;

    if (currentScore <= scoreAtRoundStart) {
      cout << "No improvement in round " << roundNum << ", stopping.\n\n" << flush;
      break;
    } else {
      cout << "Improved by " << (currentScore - scoreAtRoundStart) << ", starting another round.\n\n" << flush;
    }
  }

  result.bestScore = currentScore;
  result.bestState = currentState;

  // Update best_states file: overwrite if better, skip if worse, merge if same score. States stored alphabetically.
  // Each line: state_digits path_length state_0 state_1 ... state_n. Reverse BFS fills path when missing.
  {
    if (bestCanonicalStates.empty())
      bestCanonicalStates.insert(canonical(result.bestState));

    if (result.bestScore < fileBestScore) {
      // Worse than file: do nothing
    } else {
      ankerl::unordered_dense::map<string, BestStateEntry> byState;
      if (result.bestScore > fileBestScore) {
        for (StateType s : bestCanonicalStates)
          byState[stateToLine(s)] = BestStateEntry{stateToLine(s), -1, {}};
      } else {
        for (const BestStateEntry &e : fileEntries)
          byState[e.stateStr] = e;
        for (StateType s : bestCanonicalStates) {
          string line = stateToLine(s);
          if (byState.find(line) == byState.end())
            byState[line] = BestStateEntry{line, -1, {}};
        }
      }
      vector<BestStateEntry> toWrite;
      toWrite.reserve(byState.size());
      for (auto &kv : byState)
        toWrite.push_back(std::move(kv.second));
      sort(toWrite.begin(), toWrite.end(),
           [](const BestStateEntry &a, const BestStateEntry &b) { return a.stateStr < b.stateStr; });
      if ((int)toWrite.size() > SLIDING_WINDOW_MAX_BEST_STATES)
        toWrite.resize(SLIDING_WINDOW_MAX_BEST_STATES);

      StateType pathStartState = PATH_FROM_EMPTY ? StateType(0) : g_baseState;
      int needPath = 0;
      for (const BestStateEntry &e : toWrite)
        if (e.pathLen < 0 || RECOMPUTE_BEST_STATE_PATHS)
          needPath++;
      if (needPath > 0)
        cout << "Computing paths for " << needPath << " best state(s) via reverse BFS...\n" << flush;
      int missingPaths = 0;
      int pathIndex = 0;
      for (BestStateEntry &e : toWrite) {
        if (e.pathLen >= 0 && !RECOMPUTE_BEST_STATE_PATHS)
          continue;
        pathIndex++;
        if (needPath > 1)
          cout << "  Best-state path " << pathIndex << "/" << needPath << ": running reverse BFS...\n" << flush;
        StateType goal = lineToState(e.stateStr);
        vector<pair<int, int>> buildPath = buildPathFromReverseBfs(goal);
        if (buildPath.empty())
          continue;
        e.pathLen = (int)buildPath.size();
        vector<StateType> seq = stateSequenceFromSteps(pathStartState, buildPath);
        e.stateSequence.clear();
        e.stateSequence.reserve(seq.size());
        for (StateType st : seq)
          e.stateSequence.push_back(stateToLine(st));
        missingPaths++;
      }
      if (missingPaths > 0)
        cout << "Computed paths for " << missingPaths << " state(s) via reverse BFS.\n" << flush;

      ofstream bestStatesFile(bestStatesFilename, ios::binary);
      if (!bestStatesFile) {
        cerr << "Failed to open best-states file: " << bestStatesFilename << "\n";
      } else {
        bestStatesFile << "score " << result.bestScore << "\n";
        bestStatesFile << "count " << toWrite.size() << "\n";
        for (const BestStateEntry &e : toWrite) {
          if (e.pathLen >= 0) {
            bestStatesFile << e.stateStr << " " << e.pathLen;
            for (const string &st : e.stateSequence)
              bestStatesFile << " " << st;
            bestStatesFile << "\n";
          } else {
            bestStatesFile << e.stateStr << "\n";
          }
        }
        bestStatesFile.close();
        cout << "Wrote " << toWrite.size() << " best-score state(s) to " << bestStatesFilename << "\n" << flush;
      }
    }
  }

  cout << "\nFinal score after " << roundNum << " round(s): " << result.bestScore << "\n" << flush;

  if (SEARCH_MODE != 3) {
    // ── Build path reconstruction: reverse BFS ─────────────────────────────
    printPhaseHeader("Path reconstruction (reverse BFS)");
    printOperation("Building path from goal to start", "");
    auto pathBuildStart = chrono::steady_clock::now();
    StateType finalState = canonical(result.bestState);
    vector<pair<int, int>> buildPath = buildPathFromReverseBfs(finalState);
    auto pathBuildEnd = chrono::steady_clock::now();
    printPhaseTiming("Path reconstruction",
                     (long long)chrono::duration_cast<chrono::milliseconds>(pathBuildEnd - pathBuildStart).count());
    if (buildPath.empty()) {
      cout << "\nBuild path not found by reverse BFS.\n";
      cout << "\nFinal score: " << result.bestScore << "\n";
      result.uniqueStates = 0;
      return result;
    }

    ostream &pathOut = (g_baselineFile && g_baselineFile->is_open()) ? *g_baselineFile : cout;
    if (PATH_FROM_EMPTY)
      printExpandedPath(buildPath, StateType(0), pathOut);
    else
      printExpandedPath(expandPathForEmptyBoard(buildPath), g_baseState, pathOut);

    pathOut << "\nFinal score: " << result.bestScore << "\n";
  }
  result.uniqueStates = 0;

  return result;
}

// BFS over all reachable states.  S = uint64_t for levels 1+2, StateType for
// levels 3+.  Populates result and writes path reconstruction to cout (which
// may be redirected to the baseline file by the caller).
template <typename S> static void runBFSSearch(SearchResult &result) {
  ankerl::unordered_dense::segmented_set<S, StorageHash<S>, StorageEqual<S>> visited;

  struct Entry {
    StateType state;
    int score;
  };
  vector<Entry> q;
  size_t head = 0;

  StateType startState = PATH_FROM_EMPTY ? StateType(0) : g_baseState;
  int startScore = scoreGrid(startState);
  visited.emplace(encodeForStorage<S>(startState, {})); // sentinel: no parent
  q.push_back({startState, startScore});

  result.bestState = startState;
  result.bestScore = startScore;
  result.explored = 0;
  ankerl::unordered_dense::set<StateType, StateHash, StateEqual> bestCanonicalStates;

  while (head < q.size()) {
    auto [cur, curScore] = q[head++];
    result.explored++;

    if (curScore > result.bestScore) {
      result.bestScore = curScore;
      result.bestState = cur;
      bestCanonicalStates.clear();
      bestCanonicalStates.insert(canonical(cur));
      cout << "New best score: " << result.bestScore << "\n";
      printGrid(cur);
    } else if (curScore == result.bestScore) {
      bestCanonicalStates.insert(canonical(cur));
    }

    for (int pos = 0; pos < g_numPos; pos++) {
      int curType = stateGet(cur, pos);
      uint32_t present = neighborPresent(cur, pos);
      for (int newType = 0; newType <= g_maxType[pos]; newType++) {
        if (!FORWARD_TRANSITIONS[curType][newType])
          continue;
        if (!canPlace(present, newType))
          continue;

        auto [next, childPos] = canonicalAndChildPos(stateSet(cur, pos, newType), pos);
        int nextScore = curScore + (newType - curType);
        if (visited.emplace(encodeForStorage<S>(next, Transition{(int8_t)childPos, (int8_t)curType})).second)
          q.push_back({next, nextScore});
      }
    }

    if (result.explored % 2000000 == 0)
      cout << "Explored: " << result.explored << " | Visited: " << visited.size() << " | Best: " << result.bestScore
           << "\r" << flush;
  }

  result.uniqueStates = visited.size();

  // Write all states at best score (collected during BFS) to best_states file (with paths from BFS).
  {
    string prefix = "";
    if (std::filesystem::current_path().filename() == "build")
      prefix = "../";
    string bestStatesFilename = prefix + "best_states/level" + to_string(getCurrentLevel()) + "_type" +
                                to_string(getCurrentUnlockedType()) + ".txt";
    std::filesystem::create_directories(std::filesystem::path(bestStatesFilename).parent_path());
    StateType pathStartState = PATH_FROM_EMPTY ? StateType(0) : g_baseState;
    vector<BestStateEntry> toWrite;
    toWrite.reserve(bestCanonicalStates.size());
    for (StateType s : bestCanonicalStates) {
      vector<pair<int, int>> steps = getPathFromVisited<S>(s, visited);
      BestStateEntry e;
      e.stateStr = stateToLine(s);
      e.pathLen = (int)steps.size();
      vector<StateType> seq = stateSequenceFromSteps(pathStartState, steps);
      e.stateSequence.reserve(seq.size());
      for (StateType st : seq)
        e.stateSequence.push_back(stateToLine(st));
      toWrite.push_back(std::move(e));
    }
    sort(toWrite.begin(), toWrite.end(),
         [](const BestStateEntry &a, const BestStateEntry &b) { return a.stateStr < b.stateStr; });
    if ((int)toWrite.size() > SLIDING_WINDOW_MAX_BEST_STATES)
      toWrite.resize(SLIDING_WINDOW_MAX_BEST_STATES);
    ofstream bestStatesFile(bestStatesFilename, ios::binary);
    if (!bestStatesFile) {
      cerr << "Failed to open best-states file: " << bestStatesFilename << "\n";
    } else {
      bestStatesFile << "score " << result.bestScore << "\n";
      bestStatesFile << "count " << toWrite.size() << "\n";
      for (const BestStateEntry &e : toWrite) {
        bestStatesFile << e.stateStr << " " << e.pathLen;
        for (const string &st : e.stateSequence)
          bestStatesFile << " " << st;
        bestStatesFile << "\n";
      }
      bestStatesFile.close();
      cout << "Wrote " << toWrite.size() << " best-score state(s) to " << bestStatesFilename << "\n" << flush;
    }
  }

  ostream &pathOut = (g_baselineFile && g_baselineFile->is_open()) ? *g_baselineFile : cout;
  printBuildPath<S>(result.bestState, visited, pathOut);
}

// ── Public entrypoint for a single level/unlocked-type configuration ────────

void runSingleLevel() {
  auto totalStart = chrono::steady_clock::now();

  // (run header + unlocked type legend intentionally not printed; keep output
  // compact)
  auto setupStart = chrono::steady_clock::now();
  buildLookupTables();
  buildSymmetries();
  auto setupEnd = chrono::steady_clock::now();
  (void)setupStart;
  (void)setupEnd;

  SearchResult result{g_baseState, 0, 0, 0};
  auto phaseStart = chrono::steady_clock::now();

  if (SEARCH_MODE == 1 || SEARCH_MODE == 3) {
    result = runSlidingWindow();
  } else if (SEARCH_MODE == 2) {
    runDFSSearch(result);
  } else {
    if (usesStorage64())
      runBFSSearch<uint64_t>(result);
    else
      runBFSSearch<StateType>(result);
  }

  auto phaseEnd = chrono::steady_clock::now();
  long long phaseMs = (long long)chrono::duration_cast<chrono::milliseconds>(phaseEnd - phaseStart).count();
  const char *phaseName =
      (SEARCH_MODE == 1 || SEARCH_MODE == 3) ? "Sliding window" : (SEARCH_MODE == 2 ? "DFS" : "BFS");
  (void)phaseName;
  (void)phaseMs;

  cout << (SEARCH_MODE == 1   ? "\nSliding window complete.\n"
           : SEARCH_MODE == 2 ? "\nDFS complete.\n"
           : SEARCH_MODE == 3 ? "\nSliding window only complete (no path reconstruction).\n"
                              : "\nBFS complete.\n");
  cout << "States explored : " << result.explored << "\n";
  cout << "Unique states   : " << result.uniqueStates << "\n";
  cout << "Best score      : " << result.bestScore << "\n\n";
  cout << "Optimal layout:\n" << flush;
  printGrid(result.bestState);
  cout << flush;

  auto totalEnd = chrono::steady_clock::now();
  long long totalMs = (long long)chrono::duration_cast<chrono::milliseconds>(totalEnd - totalStart).count();
  cout << "\n";
  cout << "Total time: " << formatTimeMs(totalMs) << "\n" << flush;

  // Mode 0: path already printed by runBFSSearch. Mode 1: by runSlidingWindow.
  // Mode 2: path already printed by runDFSSearch. Mode 3: no path (sliding window only).
}
