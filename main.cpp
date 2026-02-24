#include <ankerl/unordered_dense.h>
#include <chrono>
#include <deque>
#include <functional>
#include <iostream>
#include <robin_hood.h>
#include <vector>

#include "grid.h"

using namespace std;

// Hash specialization for uint128_simd (must be after including grid.h)
namespace robin_hood {
template <> struct hash<uint128_simd> {
  size_t operator()(const uint128_simd &x) const noexcept {
    uint64_t lo = (uint64_t)x;
    uint64_t hi = (uint64_t)(x >> 64);
    return size_t(lo ^ hi);
  }
};
} // namespace robin_hood

// ── Hot-path lookup tables ─────────────────────────────────────────────────
// Flat fixed-size arrays eliminate vector pointer indirection in the inner
// loop. Max grid 7×7 = 49 positions, max 4 orthogonal neighbors each. Neighbor
// bit-shifts are precomputed (neighborIndex * 3) to avoid multiplies in
// canPlace. All positions are padded to exactly 4 entries; unused entries use
// shift 63 which reads bits 63-65 of the uint64_t — always 0 (type 0), so they
// never interfere with has[1..4] checks.

static int g_numPos;
static int g_neighborShifts[49][4]; // precomputed: neighborIndex * 3
static int g_neighborCount[49];
static int g_maxType[49];
static const uint32_t g_requiredMask[5] = {0, 1, 3, 7, 15}; // (1<<type)-1

// ── Symmetry ───────────────────────────────────────────────────────────────
// Auto-detects valid symmetries by trying all 8 dihedral transforms against
// actual grid positions. Keeps only those where every position maps to another
// valid position. Identity is always first (index 0).

static int g_numSymm;
static int g_symmPerm[8][49]; // g_symmPerm[s][i] = image of position i under
                              // symmetry s

inline StateType applyPerm(StateType s, const int *perm) {
  StateType result = 0;
  for (int i = 0; i < g_numPos; i++)
    result |= StateType((s >> (i * 3)) & 7) << (perm[i] * 3);
  return result;
}

// Canonical form = minimum StateType value over all symmetry images.
// The minimum is a unique representative of the equivalence class.
inline StateType canonical(StateType s) {
  if (g_numSymm == 1)
    return s; // No symmetries, skip computation
  StateType best = s;
  for (int k = 1; k < g_numSymm; k++) {
    StateType t = applyPerm(s, g_symmPerm[k]);
    if (t < best) {
      best = t;
      if (best == StateType(0))
        return StateType(0); // Early exit: can't get smaller than 0
    }
  }
  return best;
}

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
  cout << "Symmetries detected: " << g_numSymm << "  (search space reduced ~"
       << g_numSymm << "x)\n";
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
}

// ── State packing ──────────────────────────────────────────────────────────
// Each position stores a value 0-4 (fits in 3 bits).
// Up to 21 positions × 3 bits = 63 bits → fits in uint64_t.

inline StateType stateSet(StateType s, int pos, int val) {
  s &= ~(StateType(7) << (pos * 3));
  s |= StateType(val) << (pos * 3);
  return s;
}
inline int stateGet(StateType s, int pos) { return (s >> (pos * 3)) & 7; }

// ── Placement rule ─────────────────────────────────────────────────────────
// Tower of type T at pos requires one neighbor of each type 0..T-1.
// neighborPresent() computes the bitmask once per (state, pos); canPlace()
// then checks it cheaply for each candidate type.

inline uint32_t neighborPresent(StateType state, int pos) {
  uint32_t present = 0;
  const int *shifts = g_neighborShifts[pos];
  int count = g_neighborCount[pos];
  for (int i = 0; i < count; i++)
    present |= 1u << ((state >> shifts[i]) & 7);
  return present;
}

inline bool canPlace(uint32_t present, int type) {
  if (type == 0)
    return true;
  return (present & g_requiredMask[type]) == g_requiredMask[type];
}

// ── Path reconstruction ────────────────────────────────────────────────────
// Canonical-aware: all intermediate states are canonicalized so they match
// the transitions stored in the visited map (which live in canonical space).

template <typename MapType>
inline void printBuildPath(StateType bestState, const MapType &visited) {
  vector<pair<int, int>> steps; // (pos, newType)
  StateType cur = bestState;
  while (visited.at(cur).getPos() != -1) {
    const auto &tr = visited.at(cur);
    steps.push_back({tr.getPos(), tr.getType()});
    cur = tr.parent;
  }
  reverse(steps.begin(), steps.end());

  cout << "\n--- Build path (" << steps.size() << " steps) ---\n\n";

  StateType state = 0;
  cout << "Start:\n";
  printGrid(state);

  for (int i = 0; i < (int)steps.size(); i++) {
    int pos = steps[i].first;
    int newType = steps[i].second;
    int oldType = stateGet(state, pos);
    // Canonicalize so intermediate state matches the next transition's parent
    state = canonical(stateSet(state, pos, newType));

    cout << "\nStep " << (i + 1) << ": pos " << pos << "  " << COLORS[oldType]
         << BLACK_TEXT << COLOR_NAMES[oldType] << RESET << " (" << oldType
         << ") -> " << COLORS[newType] << BLACK_TEXT << COLOR_NAMES[newType]
         << RESET << " (" << newType << ")\n";
    printGrid(state);
  }
}

struct SearchResult {
  StateType bestState;
  int bestScore;
  long long explored;
  size_t uniqueStates;
};

// BFS from state 0 stopping as soon as target is reached.
// Returns a parent-pointer map containing the shortest path to target.
inline robin_hood::unordered_map<StateType, Transition>
bfsShortestPath(StateType target) {
  robin_hood::unordered_map<StateType, Transition> path;
  struct Entry {
    StateType state;
    int score;
  };
  vector<Entry> q;
  size_t head = 0;

  path[0] = {0, -1, -1};
  q.push_back({0, 0});

  while (head < q.size()) {
    auto [cur, curScore] = q[head++];
    if (cur == target)
      break;
    for (int pos = 0; pos < g_numPos; pos++) {
      int curType = stateGet(cur, pos);
      if (curType == 4)
        continue; // Type 4 is never a requirement, no point downgrading
      uint32_t present = neighborPresent(cur, pos);
      for (int newType = 0; newType <= g_maxType[pos]; newType++) {
        if (newType == curType)
          continue;
        if (!canPlace(present, newType))
          continue;
        StateType next = canonical(stateSet(cur, pos, newType));
        auto [it, inserted] = path.try_emplace(
            next, Transition{cur, (int8_t)pos, (int8_t)newType});
        if (inserted)
          q.push_back({next, curScore + (newType - curType)});
      }
    }
  }
  return path;
}

// Reverse BFS from a final state to base (all zeros).
// A tower can be removed if its placement requirements are met in the
// current state. The prior type at that position can be any other value,
// without needing to satisfy its own placement requirements.
static vector<pair<int, int>> buildPathFromReverseBfs(StateType finalState) {
  ankerl::unordered_dense::map<StateType, Transition> parent;
  parent.reserve(1ull << 26); // ~64M elements

  std::deque<StateType> q;

  parent[finalState] = {finalState, -1, -1};
  q.push_back(finalState);

  size_t depth = 0;
  double prevMs = 0.0;
  while (!q.empty()) {
    size_t layerCount = q.size();
    auto layerStart = chrono::steady_clock::now();
    int minTowersLeft = g_numPos;
    int maxTowersLeft = 0;
    bool anyProcessed = false;

    for (size_t i = 0; i < layerCount; i++) {
      StateType cur = q.front();
      q.pop_front();
      int towersLeft = 0;
      for (int pos = 0; pos < g_numPos; pos++) {
        if (stateGet(cur, pos) != 0)
          towersLeft++;
      }
      if (!anyProcessed) {
        minTowersLeft = towersLeft;
        maxTowersLeft = towersLeft;
        anyProcessed = true;
      } else {
        minTowersLeft = min(minTowersLeft, towersLeft);
        maxTowersLeft = max(maxTowersLeft, towersLeft);
      }
      if (cur == 0) {
        q.clear();
        break;
      }

      for (int pos = 0; pos < g_numPos; pos++) {
        int curType = stateGet(cur, pos);
        if (curType == 0)
          continue; // Nothing to remove at this position

        uint32_t present = neighborPresent(cur, pos);
        if (!canPlace(present, curType))
          continue; // This tower could not have been placed here

        for (int prevType = 0; prevType <= g_maxType[pos]; prevType++) {
          if (prevType == curType)
            continue;
          StateType next = stateSet(cur, pos, prevType);
          if (parent
                  .emplace(next, Transition{cur, (int8_t)pos, (int8_t)curType})
                  .second)
            q.push_back(next);
        }
      }
    }

    auto layerEnd = chrono::steady_clock::now();
    double layerMs =
        chrono::duration_cast<chrono::milliseconds>(layerEnd - layerStart)
            .count();
    double pctDelta =
        prevMs > 0.0 ? ((layerMs - prevMs) / prevMs) * 100.0 : 0.0;
    cout << "Reverse BFS layer " << depth << ": " << layerCount << " states, "
         << layerMs << " ms";
    if (anyProcessed) {
      cout << ", towers left min/max " << minTowersLeft << "/" << maxTowersLeft;
    }
    if (depth > 0) {
      cout << " (" << (pctDelta >= 0.0 ? "+" : "") << pctDelta << "% vs prev)";
    }
    cout << "\n";

    prevMs = layerMs;
    depth++;
  }

  if (!parent.count(0))
    return {};

  // Reconstruct forward build path: start at base state, move toward final.
  cout << "\nReverse BFS visited " << parent.size() << " states.\n";
  vector<pair<int, int>> steps;
  StateType cur = 0;
  while (parent.at(cur).getPos() != -1) {
    const auto &tr = parent.at(cur);
    steps.push_back({tr.getPos(), tr.getType()});
    cur = tr.parent;
  }
  return steps;
}

// ── Sliding Window Optimization ────────────────────────────────────────────
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

static int scoreGrid(StateType state) {
  int score = 0;
  for (int i = 0; i < g_numPos; i++) {
    int type = stateGet(state, i);
    score += type;
  }
  return score;
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
  robin_hood::unordered_set<StateType> localVisited;
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
      cout << "New best: " << bestScore << "\n";
      printGrid(bestState);
    }

    // Try modifying each position in the window
    for (int pos : windowPositions) {
      if (pos >= g_numPos)
        continue;
      int curType = stateGet(cur, pos);

      // Try all valid types for this position
      uint32_t present = neighborPresent(cur, pos);
      for (int newType = 0; newType <= g_maxType[pos]; newType++) {
        if (newType == curType)
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
  SearchResult result{0, 0, 0, 0};
  result.bestState = 0;
  result.bestScore = 0;
  result.explored = 0;

  cout << "Sliding window - sequential greedy approach:\n";
  for (const auto &row : WINDOW_KERNEL_PATTERN) {
    cout << "  " << row << "\n";
  }
  const int passOrder[] = {1, 0, 2, 3};
  auto formatPassOrder = [&]() {
    string order;
    for (size_t i = 0; i < std::size(passOrder); i++) {
      if (i > 0)
        order += "->";
      order += char('A' + passOrder[i]);
    }
    return order;
  };
  cout << "Chained 4-corner passes (" << formatPassOrder()
       << "), repeating until no improvement\n\n";

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

  StateType currentState = 0;
  int currentScore = 0;
  int roundNum = 0;

  while (true) {
    roundNum++;
    int scoreAtRoundStart = currentScore;

    cout << "============ Round " << roundNum
         << " (starting score: " << currentScore << ") ============\n\n";

    for (int passIndex = 0; passIndex < 4; passIndex++) {
      int corner = passOrder[passIndex];
      vector<int> zigzag = buildZigzagPath(corner);

      int bestScore = currentScore;
      StateType bestState = currentState;
      long long explored = 0;

      cout << "=== Pass " << (char)('A' + corner) << ": " << cornerNames[corner]
           << " (path length: " << zigzag.size() << ") ===\n";

      auto passStart = chrono::steady_clock::now();

      vector<WindowStats> stats;

      for (int windowIdx = 0; windowIdx < (int)zigzag.size(); windowIdx++) {
        vector<StateType> exploredStates;
        long long exploredBefore = explored;

        int posIdx = zigzag[windowIdx];
        auto [cr, cc] = GRID.positions[posIdx];

        auto t0 = chrono::steady_clock::now();

        exploreWindow(currentState, posIdx, exploredStates, bestScore,
                      bestState, explored);

        auto t1 = chrono::steady_clock::now();
        auto ms = chrono::duration_cast<chrono::milliseconds>(t1 - t0).count();

        bool improved = (bestState != currentState);
        if (improved) {
          currentState = bestState;
        }

        stats.push_back({posIdx, cr, cc, (int)exploredStates.size(),
                         explored - exploredBefore, ms, improved, bestScore});
      }

      currentScore = bestScore;

      auto passEnd = chrono::steady_clock::now();
      auto passMs =
          chrono::duration_cast<chrono::milliseconds>(passEnd - passStart)
              .count();

      // Print summary table
      cout << "\n┌─────┬─────────┬──────────────┬────────────────┬──────────┬──"
              "───┬───────┐\n";
      cout << "│ Pos │  (r,c)  │    Kernel    │     States     │   Time   │ "
              "Imp │ Score │\n";
      cout << "├─────┼─────────┼──────────────┼────────────────┼──────────┼────"
              "─┼───────┤\n";
      for (const auto &s : stats) {
        printf("│ %3d │ (%d, %d)  │ %12d │ %14lld │ %6lldms │  %s  │ %5d │\n",
               s.posIdx, s.row, s.col, s.kernelSize, s.statesExplored, s.timeMs,
               s.improved ? "Y" : "N", s.scoreAfter);
      }
      cout << "└─────┴─────────┴──────────────┴────────────────┴──────────┴────"
              "─┴───────┘\n";
      cout << "Pass " << (char)('A' + corner) << " (" << cornerNames[corner]
           << "): score " << currentScore << " (" << passMs << "ms, "
           << explored << " states)\n\n";

      result.explored += explored;
    }

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

  cout << "Final score after " << roundNum << " round(s): " << result.bestScore
       << "\n";

  // ── Build path reconstruction: reverse BFS ─────────────────────────────
  StateType finalState = result.bestState;
  vector<pair<int, int>> buildPath = buildPathFromReverseBfs(finalState);
  if (buildPath.empty()) {
    cout << "\nBuild path not found by reverse BFS.\n";
    cout << "\nFinal score: " << result.bestScore << "\n";
    result.uniqueStates = 0;
    return result;
  }

  cout << "\n--- Build path (" << buildPath.size() << " steps) ---\n\n";
  StateType state = 0;
  int validSteps = 0, invalidSteps = 0;
  cout << "Start:\n";
  printGrid(state);
  for (int i = 0; i < (int)buildPath.size(); i++) {
    int pos = buildPath[i].first;
    int newType = buildPath[i].second;
    int oldType = stateGet(state, pos);

    uint32_t present = neighborPresent(state, pos);
    bool valid = canPlace(present, newType);
    if (valid)
      validSteps++;
    else
      invalidSteps++;

    state = stateSet(state, pos, newType);

    auto [r, c] = GRID.positions[pos];
    cout << "\nStep " << (i + 1) << ": pos " << pos << " (" << r << "," << c
         << ")  " << COLORS[oldType] << BLACK_TEXT << COLOR_NAMES[oldType]
         << RESET << " (" << oldType << ") -> " << COLORS[newType] << BLACK_TEXT
         << COLOR_NAMES[newType] << RESET << " (" << newType << ")";
    if (!valid) {
      cout << "  \033[31m✗ INVALID\033[0m (need 0x" << hex
           << g_requiredMask[newType] << ", have 0x" << present << dec << ")";
    } else {
      cout << "  \033[32m✓\033[0m";
    }
    cout << "\n";
    printGrid(state);
  }

  int reachedScore = scoreGrid(state);
  cout << "\nBuild path: " << buildPath.size() << " steps, " << validSteps
       << " valid, " << invalidSteps << " invalid";
  if (invalidSteps == 0) {
    cout << " \033[32m— ALL VALID\033[0m";
  } else {
    cout << " \033[31m— HAS INVALID STEPS\033[0m";
  }
  cout << "\nReached score: " << reachedScore
       << " (target: " << result.bestScore << ")\n";

  cout << "\nFinal score: " << result.bestScore << "\n";
  result.uniqueStates = 0;
  return result;
}

static SearchResult runIDDFS() {
  SearchResult result{0, 0, 0, 1};
  robin_hood::unordered_set<StateType> globalSeen;
  globalSeen.insert(0);

  for (int depthLimit = 1;; depthLimit++) {
    robin_hood::unordered_set<StateType> seen;
    struct Entry {
      StateType state;
      int score;
      int depth;
    };
    vector<Entry> stack;
    size_t newStates = 0;

    stack.push_back({0, 0, 0});
    seen.insert(0);

    while (!stack.empty()) {
      auto cur = stack.back();
      stack.pop_back();
      result.explored++;

      if (cur.score > result.bestScore) {
        result.bestScore = cur.score;
        result.bestState = cur.state;
        cout << "New best score: " << result.bestScore << "\n";
        printGrid(cur.state);
      }

      if (cur.depth >= depthLimit)
        continue;

      for (int pos = 0; pos < g_numPos; pos++) {
        int curType = stateGet(cur.state, pos);
        if (curType == 4)
          continue; // Type 4 is never a requirement, no point downgrading
        uint32_t present = neighborPresent(cur.state, pos);
        for (int newType = 0; newType <= g_maxType[pos]; newType++) {
          if (newType == curType)
            continue;
          if (!canPlace(present, newType))
            continue;

          StateType next = canonical(stateSet(cur.state, pos, newType));
          int nextScore = cur.score + (newType - curType);
          if (seen.insert(next).second) {
            if (globalSeen.insert(next).second)
              newStates++;
            stack.push_back({next, nextScore, cur.depth + 1});
          }
        }
      }
    }

    result.uniqueStates = globalSeen.size();
    cout << "Depth " << depthLimit << " complete. New states: " << newStates
         << "\n";
    if (newStates == 0)
      break;
  }

  return result;
}

// ── Main ───────────────────────────────────────────────────────────────────

int main() {
  cout << "Tower Blox Deluxe Solver\n";
  cout << "========================\n";
  cout << "Tower types unlocked (0-" << MAX_UNLOCKED_TYPE << "):\n";
  for (int i = 0; i <= MAX_UNLOCKED_TYPE; i++)
    cout << "  " << COLORS[i] << BLACK_TEXT << i << RESET << " ("
         << COLOR_NAMES[i] << ")\n";
  if (SEARCH_MODE == 0)
    cout << "\nRunning BFS over all reachable states...\n\n";
  else if (SEARCH_MODE == 1)
    cout << "\nRunning DFS over all reachable states...\n\n";
  else if (SEARCH_MODE == 2)
    cout << "\nRunning IDDFS over all reachable states...\n\n";
  else
    cout << "\nRunning sliding window optimization (sequential greedy)...\n\n";

  buildLookupTables();
  buildSymmetries();

  SearchResult result{0, 0, 0, 0};
  robin_hood::unordered_map<StateType, Transition> visited;

  if (SEARCH_MODE == 3) {
    result = runSlidingWindow();
  } else if (SEARCH_MODE == 2) {
    result = runIDDFS();
  } else {
    // BFS/DFS: every state is visited exactly once.
    // For BFS, the vector is used as a queue with a head index.
    // For DFS, the vector is used as a stack via push_back/pop_back.
    struct Entry {
      StateType state;
      int score;
    };
    vector<Entry> q;
    size_t head = 0;

    StateType startState = 0;
    visited[startState] = {0, -1, -1};
    q.push_back({startState, 0});

    result.bestState = startState;
    result.bestScore = 0;
    result.explored = 0;

    while ((SEARCH_MODE == 0) ? (head < q.size()) : !q.empty()) {
      Entry curEntry;
      if (SEARCH_MODE == 0) {
        curEntry = q[head++];
      } else {
        curEntry = q.back();
        q.pop_back();
      }
      auto [cur, curScore] = curEntry;
      result.explored++;

      if (curScore > result.bestScore) {
        result.bestScore = curScore;
        result.bestState = cur;
        cout << "New best score: " << result.bestScore << "\n";
        printGrid(cur);
      }

      // Try every possible single-position change
      for (int pos = 0; pos < g_numPos; pos++) {
        int curType = stateGet(cur, pos);
        if (curType == 4)
          continue; // Type 4 is never a requirement, no point downgrading
        uint32_t present = neighborPresent(cur, pos); // computed once per pos
        for (int newType = 0; newType <= g_maxType[pos]; newType++) {
          if (newType == curType)
            continue;
          if (!canPlace(present, newType))
            continue;

          StateType next = canonical(stateSet(cur, pos, newType));
          int nextScore = curScore + (newType - curType);
          auto [it, inserted] = visited.try_emplace(
              next, Transition{cur, (int8_t)pos, (int8_t)newType});
          if (inserted) {
            q.push_back({next, nextScore});
          }
        }
      }

      if (result.explored % 2000000 == 0)
        cout << "Explored: " << result.explored
             << " | Visited: " << visited.size()
             << " | Best: " << result.bestScore << "\r" << flush;
    }

    result.uniqueStates = visited.size();
  }

  if (SEARCH_MODE == 3)
    cout << "\nSliding window complete.\n";
  else if (SEARCH_MODE == 2)
    cout << "\nIDDFS complete.\n";
  else if (SEARCH_MODE == 1)
    cout << "\nDFS complete.\n";
  else
    cout << "\nBFS complete.\n";
  cout << "States explored : " << result.explored << "\n";
  cout << "Unique states   : " << result.uniqueStates << "\n";
  cout << "Best score      : " << result.bestScore << "\n\n";
  cout << "Optimal layout:\n";
  printGrid(result.bestState);

  // Path reconstruction only available for BFS
  if (SEARCH_MODE == 0) {
    printBuildPath(result.bestState, visited);
  } else if (SEARCH_MODE == 3) {
    // Path already printed by runSlidingWindow
  } else {
    cout << "\nFinding shortest build path via BFS...\n";
    auto pathMap = bfsShortestPath(result.bestState);
    printBuildPath(result.bestState, pathMap);
  }

  return 0;
}
