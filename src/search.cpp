#include <ankerl/unordered_dense.h>
#include <chrono>
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

// Global file stream for baseline output (only used for reverse BFS section)
ofstream *g_baselineFile = nullptr;
streambuf *g_originalCoutBuf = nullptr;

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
static int scoreGrid(StateType state) {
  int score = 0;
  for (int i = 0; i < g_numPos; i++)
    score += stateGet(state, i) - 1;
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
static vector<pair<int, int>>
expandPathForEmptyBoard(const vector<pair<int, int>> &steps) {
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

static void advanceStateForDisplay(StateType &canonState,
                                   StateType &displayState,
                                   uint64_t &placedMask, int pos, int newType) {
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

static void printExpandedPath(const vector<pair<int, int>> &expandedSteps,
                              StateType startState) {
  const bool fromEmpty = (startState == StateType(0));
  StateType canonState = startState;
  StateType displayState = startState;
  uint64_t placedMask = 0;

  // Pre-pass to compute total steps (and trailing implicit blues when from
  // all-1).
  for (auto [pos, newType] : expandedSteps)
    advanceStateForDisplay(canonState, displayState, placedMask, pos, newType);

  int totalSteps = (int)expandedSteps.size();
  if (!fromEmpty) {
    int missingCount = 0;
    for (int i = 0; i < g_numPos; i++) {
      if (!((placedMask >> i) & 1))
        missingCount++;
    }
    totalSteps += missingCount;
  }

  // Reset for actual print pass.
  canonState = startState;
  displayState = startState;
  placedMask = 0;

  cout << "\n--- Build path (" << totalSteps << " steps) ---\n\n";
  cout << "Start (" << (fromEmpty ? "empty" : "all blue") << " board):\n";
  printGrid(startState, uint64_t(0));
  int validSteps = 0, invalidSteps = 0;
  int stepNum = 1;

  for (auto [pos, newType] : expandedSteps) {
    // pos is a position index in canonical space.
    // Find which symmetry maps our displayState to its canonical form,
    // then invert it to find the physical position.
    int usedSymm = 0;
    for (int k = 0; k < g_numSymm; k++) {
      if (applyPerm(displayState, g_symmPerm[k]) == canonState) {
        usedSymm = k;
        break;
      }
    }
    // Build inverse permutation: inversePerm[canonical_pos] = physical_pos
    int inversePerm[49];
    for (int i = 0; i < g_numPos; i++)
      inversePerm[g_symmPerm[usedSymm][i]] = i;

    int physPos = inversePerm[pos];
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
    cout << "\nStep " << stepNum << ": pos " << physPos << " (" << r << "," << c
         << ")  ";
    if (oldType == -1) {
      cout << "empty -> " << COLORS[newType] << BLACK_TEXT
           << COLOR_NAMES[newType] << RESET << " (" << newType << ")";
    } else {
      cout << COLORS[oldType] << BLACK_TEXT << COLOR_NAMES[oldType] << RESET
           << " (" << oldType << ") -> " << COLORS[newType] << BLACK_TEXT
           << COLOR_NAMES[newType] << RESET << " (" << newType << ")";
    }
    cout << (valid ? "  \033[32m\u2713\033[0m"
                   : "  \033[31m\u2717 INVALID\033[0m")
         << "\n";
    printGrid(displayState, placedMask);
    stepNum++;
  }

  // When path started from all-blue: any position never explicitly placed is an
  // implicit blue; show a final step for each.
  if (!fromEmpty) {
    for (int i = 0; i < g_numPos; i++) {
      if (!((placedMask >> i) & 1)) {
        placedMask |= (1ull << i);
        auto [r, c] = GRID.positions[i];
        cout << "\nStep " << stepNum << ": pos " << i << " (" << r << "," << c
             << ")  ";
        cout << "empty -> " << COLORS[1] << BLACK_TEXT << COLOR_NAMES[1]
             << RESET << " (1)";
        cout << "  \033[32m\u2713\033[0m"
             << "\n";
        printGrid(displayState, placedMask);
        stepNum++;
        validSteps++;
      }
    }
  }

  int reachedScore = scoreGrid(displayState);
  cout << "\nBuild path: " << totalSteps << " steps, " << validSteps
       << " valid, " << invalidSteps << " invalid";
  if (invalidSteps == 0) {
    cout << " \033[32m\u2014 ALL VALID\033[0m";
  } else {
    cout << " \033[31m\u2014 HAS INVALID STEPS\033[0m";
  }
  cout << "\nReached score: " << reachedScore << "\n";
}

template <typename S>
inline void printBuildPath(StateType bestState,
                           const ankerl::unordered_dense::segmented_set<
                               S, StorageHash<S>, StorageEqual<S>> &visited) {
  vector<pair<int, int>>
      steps; // (pos, newType) — pos in state *before* the move
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
    printExpandedPath(steps, StateType(0));
  else
    printExpandedPath(expandPathForEmptyBoard(steps), g_baseState);
}

void printBuildPathSteps(const vector<pair<int, int>> &buildPath) {
  if (PATH_FROM_EMPTY)
    printExpandedPath(buildPath, StateType(0));
  else
    printExpandedPath(expandPathForEmptyBoard(buildPath), g_baseState);
}

// BFS from base state (all blue) stopping as soon as target is reached.
// Returns a set of encoded states (transition packed into top byte) for path
// reconstruction.
inline ankerl::unordered_dense::segmented_set<StateType, StateHash, StateEqual>
bfsShortestPath(StateType target) {
  ankerl::unordered_dense::segmented_set<StateType, StateHash, StateEqual> path;
  struct Entry {
    StateType state;
    int score;
  };
  vector<Entry> q;
  size_t head = 0;

  StateType bfsStart = PATH_FROM_EMPTY ? StateType(0) : g_baseState;
  path.emplace(encodeTransition(bfsStart, {})); // sentinel: packed=255
  q.push_back({bfsStart, 0});

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
        auto [next, childPos] =
            canonicalAndChildPos(stateSet(cur, pos, newType), pos);
        if (path.emplace(encodeTransition(next, Transition{(int8_t)childPos,
                                                           (int8_t)curType}))
                .second)
          q.push_back({next, curScore + (newType - curType)});
      }
    }
  }
  return path;
}

// ── Sliding Window Optimization ─────────────────────────────────────────────-
// Processes grid using 3x3 windows with stride=1 along a JPEG-like zigzag path
// Each window respects constraints from the entire grid (in+out of window)

// Build zigzag path starting from a given corner
// corner: 0=top-left, 1=top-right, 2=bottom-left, 3=bottom-right
static vector<int> buildZigzagPath(int corner = 0) {
  int R = GRID.rows, C = GRID.cols;
  vector<vector<pair<int, int>>> diagonals;

  // Group by diagonal - direction depends on corner
  for (int r = 0; r < R; r++) {
    for (int c = 0; c < C; c++) {
      pair<int, int> coord = {r, c};
      if (!GRID.coordToIdx.count(coord))
        continue;

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
    // Alternate direction for each diagonal (zigzag effect)
    if (d % 2 == 1)
      reverse(diag.begin(), diag.end());
    for (auto [r, c] : diag) {
      path.push_back(GRID.coordToIdx[{r, c}]);
    }
  }
  return path;
}

// BFS within a window: exhaustively explore ALL reachable states
// by modifying only positions within the window
static void exploreWindow(StateType startState, int windowCenterIdx,
                          vector<StateType> &newStates, int &bestScore,
                          StateType &bestState, long long &explored) {
  // Get kernel pattern window around position
  auto [cr, cc] = GRID.positions[windowCenterIdx];

  // Parse kernel pattern to get offsets
  int patternRows = (int)WINDOW_KERNEL_PATTERN.size();
  int patternCols = 0;
  for (const auto &row : WINDOW_KERNEL_PATTERN) {
    patternCols = max(patternCols, (int)row.size());
  }
  int centerR = patternRows / 2;
  int centerC = patternCols / 2;

  vector<int> windowPositions;
  for (int pr = 0; pr < patternRows; pr++) {
    for (int pc = 0; pc < (int)WINDOW_KERNEL_PATTERN[pr].size(); pc++) {
      if (WINDOW_KERNEL_PATTERN[pr][pc] == '.') {
        int dr = pr - centerR;
        int dc = pc - centerC;
        auto it = GRID.coordToIdx.find({cr + dr, cc + dc});
        if (it != GRID.coordToIdx.end()) {
          windowPositions.push_back(it->second);
        }
      }
    }
  }

  // Full BFS within window - explore ALL states reachable by modifying window
  // positions
  ankerl::unordered_dense::segmented_set<StateType, StateHash, StateEqual>
      localVisited;
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

    // Stop expanding this branch if we've reached max depth
    if (depth >= MAX_DEPTH)
      continue;

    if (curScore > bestScore) {
      bestScore = curScore;
      bestState = cur;
      if (g_activeStatusLine)
        g_activeStatusLine->clear();
      cout << "New best: " << bestScore << "\n";
      printGrid(bestState);
      if (g_activeStatusLine)
        g_activeStatusLine->redraw();
    }

    // Try modifying each position in the window
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

        // Only track within this window to prevent infinite loops
        if (localVisited.insert(next).second) {
          int nextScore = curScore + (newType - curType);
          q.push_back({next, nextScore, depth + 1});
          newStates.push_back(next);
        }
      }
    }
  }
  // localVisited is destroyed here - no global state persists
}

static SearchResult runSlidingWindow() {
  SearchResult result{g_baseState, 0, 0, 0};
  result.bestState = g_baseState;
  result.bestScore = 0;
  result.explored = 0;

  // (kernel pattern intentionally not printed; keep output compact)

  const char *cornerNames[] = {"top-left", "top-right", "bottom-left",
                               "bottom-right"};

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
  int currentScore = 0;
  int roundNum = 0;

  while (true) {
    roundNum++;
    int scoreAtRoundStart = currentScore;

    // (round header intentionally not printed; keep output compact)

    int corner = 1; // top-right
    vector<int> zigzag = buildZigzagPath(corner);

    int bestScore = currentScore;
    StateType bestState = currentState;
    long long explored = 0;

    int totalWindows = (int)zigzag.size();
    // (pass header intentionally not printed; keep output compact)

    auto passStart = chrono::steady_clock::now();

    vector<WindowStats> stats;
    DynamicStatusLine passLine(g_originalCoutBuf == nullptr);
    g_activeStatusLine = &passLine;

    for (int windowIdx = 0; windowIdx < totalWindows; windowIdx++) {
      vector<StateType> exploredStates;
      long long exploredBefore = explored;

      int posIdx = zigzag[windowIdx];
      auto [cr, cc] = GRID.positions[posIdx];

      {
        string passName = string(1, (char)('A' + corner));
        string win = padLeft(to_string(windowIdx + 1), 2) + "/" +
                     padLeft(to_string(totalWindows), 2);
        string center = "(" + padLeft(to_string(cr), 2) + "," +
                        padLeft(to_string(cc), 2) + ")";
        string line =
            "  Pass " + passName + " " +
            progressBar((size_t)(windowIdx + 1), (size_t)totalWindows) +
            " | win " + win + " | center " + center + " | best " +
            padLeft(to_string(bestScore), 3) + " | explored " +
            padLeft(formatNum(explored), 12);
        passLine.update(line);
      }
      auto t0 = chrono::steady_clock::now();

      exploreWindow(currentState, posIdx, exploredStates, bestScore, bestState,
                    explored);

      auto t1 = chrono::steady_clock::now();
      auto ms = chrono::duration_cast<chrono::milliseconds>(t1 - t0).count();

      bool improved = (bestState != currentState);
      if (improved) {
        currentState = bestState;
      }

      long long deltaStates = explored - exploredBefore;
      {
        string passName = string(1, (char)('A' + corner));
        string win = padLeft(to_string(windowIdx + 1), 2) + "/" +
                     padLeft(to_string(totalWindows), 2);
        string center = "(" + padLeft(to_string(cr), 2) + "," +
                        padLeft(to_string(cc), 2) + ")";
        string line =
            "  Pass " + passName + " " +
            progressBar((size_t)(windowIdx + 1), (size_t)totalWindows) +
            " | win " + win + " | center " + center + " | states " +
            padLeft(formatNum(deltaStates), 12) + " | step " +
            padLeft(formatTimeMs((long long)ms), 8) + " | best " +
            padLeft(to_string(bestScore), 3) + " | explored " +
            padLeft(formatNum(explored), 12);
        passLine.update(line);
      }

      stats.push_back({posIdx, cr, cc, (int)exploredStates.size(),
                       explored - exploredBefore, ms, improved, bestScore});
    }
    passLine.finish();
    g_activeStatusLine = nullptr;

    currentScore = bestScore;

    auto passEnd = chrono::steady_clock::now();
    auto passMs =
        chrono::duration_cast<chrono::milliseconds>(passEnd - passStart)
            .count();

    // Print summary table (built from widths; no hardcoded border strings)
    cout << "\n";
    const vector<int> widths = {5, 9, 14, 16, 12, 6, 7};
    const vector<TableAlign> headerAlign = {
        TableAlign::Center, TableAlign::Center, TableAlign::Center,
        TableAlign::Center, TableAlign::Center, TableAlign::Center,
        TableAlign::Center};
    const vector<TableAlign> rowAlign = {TableAlign::Right, TableAlign::Center,
                                         TableAlign::Right, TableAlign::Right,
                                         TableAlign::Right, TableAlign::Center,
                                         TableAlign::Right};

    printTableBorderTop(widths);
    printTableRow({"Pos", "(r,c)", "Kernel", "States", "Time", "Imp", "Score"},
                  widths, headerAlign);
    printTableBorderMid(widths);
    for (const auto &s : stats) {
      printTableRow({to_string(s.posIdx),
                     "(" + to_string(s.row) + ", " + to_string(s.col) + ")",
                     formatNum((size_t)s.kernelSize),
                     formatNum((long long)s.statesExplored),
                     formatTimeMs(s.timeMs), (s.improved ? "Y" : "N"),
                     to_string(s.scoreAfter)},
                    widths, rowAlign);
    }
    printTableBorderBottom(widths);
    cout << "Pass " << (char)('A' + corner) << " (" << cornerNames[corner]
         << "): score " << currentScore << " ("
         << formatTimeMs((long long)passMs) << ", " << formatNum(explored)
         << " states)\n\n";

    result.explored += explored;

    cout << "Round " << roundNum << " complete: score " << scoreAtRoundStart
         << " -> " << currentScore << "\n";

    if (currentScore <= scoreAtRoundStart) {
      cout << "No improvement in round " << roundNum << ", stopping.\n\n";
      break;
    } else {
      cout << "Improved by " << (currentScore - scoreAtRoundStart)
           << ", starting another round.\n\n";
    }
  }

  result.bestScore = currentScore;
  result.bestState = currentState;

  cout << "\nFinal score after " << roundNum
       << " round(s): " << result.bestScore << "\n";

  // ── Build path reconstruction: reverse BFS ─────────────────────────────
  printPhaseHeader("Path reconstruction (reverse BFS)");
  printOperation("Building path from goal to start", "");
  auto pathBuildStart = chrono::steady_clock::now();
  StateType finalState = canonical(result.bestState);
  vector<pair<int, int>> buildPath = buildPathFromReverseBfs(finalState);
  auto pathBuildEnd = chrono::steady_clock::now();
  printPhaseTiming("Path reconstruction",
                   (long long)chrono::duration_cast<chrono::milliseconds>(
                       pathBuildEnd - pathBuildStart)
                       .count());
  if (buildPath.empty()) {
    cout << "\nBuild path not found by reverse BFS.\n";
    cout << "\nFinal score: " << result.bestScore << "\n";
    result.uniqueStates = 0;
    return result;
  }

  if (PATH_FROM_EMPTY)
    printExpandedPath(buildPath, StateType(0));
  else
    printExpandedPath(expandPathForEmptyBoard(buildPath), g_baseState);

  cout << "\nFinal score: " << result.bestScore << "\n";
  result.uniqueStates = 0;

  // Restore cout if it was redirected
  if (g_originalCoutBuf) {
    cout.rdbuf(g_originalCoutBuf);
    g_originalCoutBuf = nullptr;
  }

  return result;
}

// BFS over all reachable states.  S = uint64_t for levels 1+2, StateType for
// levels 3+.  Populates result and writes path reconstruction to cout (which
// may be redirected to the baseline file by the caller).
template <typename S> static void runBFSSearch(SearchResult &result) {
  ankerl::unordered_dense::segmented_set<S, StorageHash<S>, StorageEqual<S>>
      visited;

  struct Entry {
    StateType state;
    int score;
  };
  vector<Entry> q;
  size_t head = 0;

  StateType startState = PATH_FROM_EMPTY ? StateType(0) : g_baseState;
  visited.emplace(encodeForStorage<S>(startState, {})); // sentinel: no parent
  q.push_back({startState, 0});

  result.bestState = startState;
  result.bestScore = 0;
  result.explored = 0;

  while (head < q.size()) {
    auto [cur, curScore] = q[head++];
    result.explored++;

    if (curScore > result.bestScore) {
      result.bestScore = curScore;
      result.bestState = cur;
      cout << "New best score: " << result.bestScore << "\n";
      printGrid(cur);
    }

    for (int pos = 0; pos < g_numPos; pos++) {
      int curType = stateGet(cur, pos);
      uint32_t present = neighborPresent(cur, pos);
      for (int newType = 0; newType <= g_maxType[pos]; newType++) {
        if (!FORWARD_TRANSITIONS[curType][newType])
          continue;
        if (!canPlace(present, newType))
          continue;

        auto [next, childPos] =
            canonicalAndChildPos(stateSet(cur, pos, newType), pos);
        int nextScore = curScore + (newType - curType);
        if (visited
                .emplace(encodeForStorage<S>(
                    next, Transition{(int8_t)childPos, (int8_t)curType}))
                .second)
          q.push_back({next, nextScore});
      }
    }

    if (result.explored % 2000000 == 0)
      cout << "Explored: " << result.explored
           << " | Visited: " << visited.size()
           << " | Best: " << result.bestScore << "\r" << flush;
  }

  result.uniqueStates = visited.size();

  if (g_baselineFile && g_baselineFile->is_open()) {
    g_originalCoutBuf = cout.rdbuf();
    cout.rdbuf(g_baselineFile->rdbuf());
  }
  printBuildPath<S>(result.bestState, visited);
  if (g_originalCoutBuf) {
    cout.rdbuf(g_originalCoutBuf);
    g_originalCoutBuf = nullptr;
  }
}

// ── Public entrypoint for a single level/unlocked-type configuration ────────

void runSingleLevel() {
  auto totalStart = chrono::steady_clock::now();

  // (run header + unlocked type legend intentionally not printed; keep output compact)
  auto setupStart = chrono::steady_clock::now();
  buildLookupTables();
  buildSymmetries();
  auto setupEnd = chrono::steady_clock::now();
  (void)setupStart;
  (void)setupEnd;

  SearchResult result{g_baseState, 0, 0, 0};
  auto phaseStart = chrono::steady_clock::now();

  if (SEARCH_MODE == 1) {
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
  long long phaseMs = (long long)chrono::duration_cast<chrono::milliseconds>(
                          phaseEnd - phaseStart)
                          .count();
  const char *phaseName =
      SEARCH_MODE == 1 ? "Sliding window" : (SEARCH_MODE == 2 ? "DFS" : "BFS");
  (void)phaseName;
  (void)phaseMs;

  cout << (SEARCH_MODE == 1   ? "\nSliding window complete.\n"
           : SEARCH_MODE == 2 ? "\nDFS complete.\n"
                              : "\nBFS complete.\n");
  cout << "States explored : " << result.explored << "\n";
  cout << "Unique states   : " << result.uniqueStates << "\n";
  cout << "Best score      : " << result.bestScore << "\n\n";
  cout << "Optimal layout:\n";
  printGrid(result.bestState);

  auto totalEnd = chrono::steady_clock::now();
  long long totalMs = (long long)chrono::duration_cast<chrono::milliseconds>(
                          totalEnd - totalStart)
                          .count();
  cout << "\n";
  cout << "Total time: " << formatTimeMs(totalMs) << "\n";

  // Mode 0: path already printed by runBFSSearch. Mode 1: by runSlidingWindow.
  // Mode 2: path already printed by runDFSSearch.
}
