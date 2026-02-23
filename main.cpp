#include <chrono>
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

// ── Sliding Window Optimization ────────────────────────────────────────────
// Processes grid using 3x3 windows with stride=1 along a JPEG-like zigzag path
// Each window respects constraints from the entire grid (in+out of window)

static vector<int> buildZigzagPath() {
  // Build diagonal/zigzag path through grid positions
  // Groups positions by diagonal, then alternates direction per diagonal
  int R = GRID.rows, C = GRID.cols;
  vector<vector<pair<int, int>>> diagonals;
  
  // Group by diagonal (r + c)
  for (int r = 0; r < R; r++) {
    for (int c = 0; c < C; c++) {
      pair<int, int> coord = {r, c};
      if (!GRID.coordToIdx.count(coord)) continue;
      int diag = r + c;
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
                             vector<StateType> &newStates,
                             int &bestScore, StateType &bestState,
                             long long &explored, bool verbose = false) {
  // Get kernel pattern window around position
  auto [cr, cc] = GRID.positions[windowCenterIdx];
  
  if (verbose) {
    cout << "  Window at position " << windowCenterIdx << " (row " << cr << ", col " << cc << ")";
  }
  
  size_t statesBeforeWindow = newStates.size();
  long long movesBefore = explored;
  
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
  
  // Debug: show kernel details
  if (verbose) {
    cout << " [kernel size=" << windowPositions.size() << ", positions={";
    for (size_t i = 0; i < windowPositions.size(); i++) {
      int pos = windowPositions[i];
      cout << pos << "(" << g_maxType[pos] << ")";
      if (i < windowPositions.size() - 1) cout << ",";
    }
    cout << "}]";
  }
  
  // Full BFS within window - explore ALL states reachable by modifying window positions
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
    if (depth >= MAX_DEPTH) continue;
    
    if (curScore > bestScore) {
      bestScore = curScore;
      bestState = cur;
      cout << "New best: " << bestScore << "\n";
      printGrid(bestState);
    }
    
    // Try modifying each position in the window
    for (int pos : windowPositions) {
      if (pos >= g_numPos) continue;
      int curType = stateGet(cur, pos);
      
      // Try all valid types for this position
      uint32_t present = neighborPresent(cur, pos);
      for (int newType = 0; newType <= g_maxType[pos]; newType++) {
        if (newType == curType) continue;
        if (!canPlace(present, newType)) continue;
        
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
  
  size_t statesGenerated = newStates.size() - statesBeforeWindow;
  size_t windowStatesExplored = localVisited.size();
  long long movesTried = explored - movesBefore;
  
  if (verbose && statesGenerated > 0) {
    cout << " -> explored " << windowStatesExplored << " states (" 
         << movesTried << " moves), generated " << statesGenerated << " new\n";
  } else if (verbose) {
    cout << " -> explored " << windowStatesExplored << " states, no new\n";
  }
}

static SearchResult runSlidingWindow() {
  SearchResult result{0, 0, 0, 0};
  vector<int> zigzag = buildZigzagPath();
  
  StateType currentState = 0;
  result.bestState = currentState;
  result.bestScore = 0;
  result.explored = 0;
  
  cout << "Sliding window - sequential greedy approach:\n";
  for (const auto &row : WINDOW_KERNEL_PATTERN) {
    cout << "  " << row << "\n";
  }
  cout << "Zigzag path length: " << zigzag.size() << "\n";
  cout << "For each position: explore all reachable states, keep best, move on\n\n";
  
  // Process each window position sequentially
  struct WindowStats {
    int posIdx;
    int row, col;
    int kernelSize;
    long long statesExplored;
    long long timeMs;
    bool improved;
    int scoreAfter;
  };
  vector<WindowStats> stats;
  
  for (int windowIdx = 0; windowIdx < (int)zigzag.size(); windowIdx++) {
    vector<StateType> exploredStates;
    long long exploredBefore = result.explored;
    
    auto [cr, cc] = GRID.positions[zigzag[windowIdx]];
    
    auto t0 = chrono::steady_clock::now();
    
    exploreWindow(currentState, windowIdx, exploredStates,
                 result.bestScore, result.bestState, result.explored, false);
    
    auto t1 = chrono::steady_clock::now();
    auto ms = chrono::duration_cast<chrono::milliseconds>(t1 - t0).count();
    
    bool improved = (result.bestState != currentState);
    if (improved) {
      currentState = result.bestState;
    }
    
    stats.push_back({zigzag[windowIdx], cr, cc, (int)exploredStates.size(),
                     result.explored - exploredBefore, ms, improved, result.bestScore});
  }
  
  // Print summary table
  cout << "\n┌─────┬─────────┬──────────────┬────────────────┬──────────┬─────┬───────┐\n";
  cout << "│ Pos │  (r,c)  │    Kernel    │     States     │   Time   │ Imp │ Score │\n";
  cout << "├─────┼─────────┼──────────────┼────────────────┼──────────┼─────┼───────┤\n";
  for (const auto &s : stats) {
    printf("│ %3d │ (%d, %d)  │ %12d │ %14lld │ %6lldms │  %s  │ %5d │\n",
           s.posIdx, s.row, s.col, s.kernelSize, s.statesExplored,
           s.timeMs, s.improved ? "Y" : "N", s.scoreAfter);
  }
  cout << "└─────┴─────────┴──────────────┴────────────────┴──────────┴─────┴───────┘\n";
  
  // Compute minimal build path from final state
  // Place towers in type order (1, 2, 3, 4) - prerequisites are always satisfied
  StateType finalState = result.bestState;
  vector<pair<int,int>> buildPath; // (pos, type)
  for (int targetType = 1; targetType <= 4; targetType++) {
    for (int pos = 0; pos < g_numPos; pos++) {
      if (stateGet(finalState, pos) == targetType) {
        buildPath.push_back({pos, targetType});
      }
    }
  }
  
  cout << "\n--- Shortest build path (" << buildPath.size() << " steps) ---\n\n";
  StateType state = 0;
  cout << "Start:\n";
  printGrid(state);
  for (int i = 0; i < (int)buildPath.size(); i++) {
    int pos = buildPath[i].first;
    int newType = buildPath[i].second;
    int oldType = stateGet(state, pos);
    state = stateSet(state, pos, newType);
    
    auto [r, c] = GRID.positions[pos];
    cout << "\nStep " << (i + 1) << ": pos " << pos << " (" << r << "," << c << ")  "
         << COLORS[oldType] << BLACK_TEXT << COLOR_NAMES[oldType] << RESET
         << " (" << oldType << ") -> "
         << COLORS[newType] << BLACK_TEXT << COLOR_NAMES[newType] << RESET
         << " (" << newType << ")\n";
    printGrid(state);
  }
  
  cout << "\nFinal score: " << result.bestScore << "\n";
  result.uniqueStates = 0; // Not tracking globally
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
