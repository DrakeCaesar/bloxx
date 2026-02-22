#include <functional>
#include <iostream>
#include <robin_hood.h>
#include <vector>

#include "grid.h"

using namespace std;

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

inline uint64_t applyPerm(uint64_t s, const int *perm) {
  uint64_t result = 0;
  for (int i = 0; i < g_numPos; i++)
    result |= uint64_t((s >> (i * 3)) & 7) << (perm[i] * 3);
  return result;
}

// Canonical form = minimum uint64_t value over all symmetry images.
// The minimum is a unique representative of the equivalence class.
inline uint64_t canonical(uint64_t s) {
  uint64_t best = s;
  for (int k = 1; k < g_numSymm; k++) {
    uint64_t t = applyPerm(s, g_symmPerm[k]);
    if (t < best)
      best = t;
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

inline uint64_t stateSet(uint64_t s, int pos, int val) {
  s &= ~(uint64_t(7) << (pos * 3));
  s |= uint64_t(val) << (pos * 3);
  return s;
}
inline int stateGet(uint64_t s, int pos) { return (s >> (pos * 3)) & 7; }

// ── Placement rule ─────────────────────────────────────────────────────────
// Tower of type T at pos requires one neighbor of each type 0..T-1.
// neighborPresent() computes the bitmask once per (state, pos); canPlace()
// then checks it cheaply for each candidate type.

inline uint32_t neighborPresent(uint64_t state, int pos) {
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
inline void printBuildPath(uint64_t bestState, const MapType &visited) {
  vector<pair<int, int>> steps; // (pos, newType)
  uint64_t cur = bestState;
  while (visited.at(cur).pos != -1) {
    const auto &tr = visited.at(cur);
    steps.push_back({tr.pos, tr.type});
    cur = tr.parent;
  }
  reverse(steps.begin(), steps.end());

  cout << "\n--- Build path (" << steps.size() << " steps) ---\n\n";

  uint64_t state = 0;
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
  uint64_t bestState;
  int bestScore;
  long long explored;
  size_t uniqueStates;
};

// BFS from state 0 stopping as soon as target is reached.
// Returns a parent-pointer map containing the shortest path to target.
inline robin_hood::unordered_map<uint64_t, Transition>
bfsShortestPath(uint64_t target) {
  robin_hood::unordered_map<uint64_t, Transition> path;
  struct Entry {
    uint64_t state;
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
        uint64_t next = canonical(stateSet(cur, pos, newType));
        auto [it, inserted] = path.try_emplace(
            next, Transition{cur, (int8_t)pos, (int8_t)newType});
        if (inserted)
          q.push_back({next, curScore + (newType - curType)});
      }
    }
  }
  return path;
}

static SearchResult runIDDFS() {
  SearchResult result{0, 0, 0, 1};
  robin_hood::unordered_set<uint64_t> globalSeen;
  globalSeen.insert(0);

  for (int depthLimit = 1;; depthLimit++) {
    robin_hood::unordered_set<uint64_t> seen;
    struct Entry {
      uint64_t state;
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

          uint64_t next = canonical(stateSet(cur.state, pos, newType));
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
  cout << "Tower types:\n";
  for (int i = 0; i <= 4; i++)
    cout << "  " << COLORS[i] << BLACK_TEXT << i << RESET << " ("
         << COLOR_NAMES[i] << ")\n";
  if (SEARCH_MODE == 0)
    cout << "\nRunning BFS over all reachable states...\n\n";
  else if (SEARCH_MODE == 1)
    cout << "\nRunning DFS over all reachable states...\n\n";
  else
    cout << "\nRunning IDDFS over all reachable states...\n\n";

  buildLookupTables();
  buildSymmetries();

  SearchResult result{0, 0, 0, 0};
  robin_hood::unordered_map<uint64_t, Transition> visited;

  if (SEARCH_MODE == 2) {
    result = runIDDFS();
  } else {
    // BFS/DFS: every state is visited exactly once.
    // For BFS, the vector is used as a queue with a head index.
    // For DFS, the vector is used as a stack via push_back/pop_back.
    struct Entry {
      uint64_t state;
      int score;
    };
    vector<Entry> q;
    size_t head = 0;

    uint64_t startState = 0;
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

          uint64_t next = canonical(stateSet(cur, pos, newType));
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

  if (SEARCH_MODE == 2)
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

  if (SEARCH_MODE == 0) {
    printBuildPath(result.bestState, visited);
  } else {
    cout << "\nFinding shortest build path via BFS...\n";
    auto pathMap = bfsShortestPath(result.bestState);
    printBuildPath(result.bestState, pathMap);
  }

  return 0;
}
