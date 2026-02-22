#include <iostream>
#include <unordered_map>
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
inline int stateScore(uint64_t s) {
  int total = 0;
  for (int i = 0; i < g_numPos; i++)
    total += (s >> (i * 3)) & 7;
  return total;
}

// ── Placement rule ─────────────────────────────────────────────────────────
// Tower of type T at pos requires one neighbor of each type 0..T-1.

inline bool canPlace(uint64_t state, int pos, int type) {
  if (type == 0)
    return true;
  // Bitmask of which types are present among neighbors (bit i = type i seen)
  uint32_t present = 0;
  const int *shifts = g_neighborShifts[pos];
  int count = g_neighborCount[pos];
  for (int i = 0; i < count; i++)
    present |= 1u << ((state >> shifts[i]) & 7);
  // Need all of types 0..type-1: use precomputed mask
  return (present & g_requiredMask[type]) == g_requiredMask[type];
}

// ── Main ───────────────────────────────────────────────────────────────────

int main() {
  cout << "Tower Blox Deluxe Solver\n";
  cout << "========================\n";
  cout << "Tower types:\n";
  for (int i = 0; i <= 4; i++)
    cout << "  " << COLORS[i] << BLACK_TEXT << i << RESET << " ("
         << COLOR_NAMES[i] << ")\n";
  cout << "\nRunning BFS over all reachable states...\n\n";

  buildLookupTables();

  // BFS: every state is visited exactly once.
  // Queue stores (state, score) pairs so score is never recomputed.
  // A flat vector with a head index is a single contiguous allocation,
  // avoiding std::deque's chunked heap blocks.
  unordered_map<uint64_t, Transition> visited;
  visited.reserve(1 << 23);

  struct Entry {
    uint64_t state;
    int score;
  };
  vector<Entry> q;
  q.reserve(1 << 23);
  size_t head = 0;

  uint64_t startState = 0;
  visited[startState] = {0, -1, -1};
  q.push_back({startState, 0});

  uint64_t bestState = startState;
  int bestScore = 0;
  long long explored = 0;

  while (head < q.size()) {
    auto [cur, curScore] = q[head++];
    explored++;

    if (curScore > bestScore) {
      bestScore = curScore;
      bestState = cur;
      cout << "New best score: " << bestScore << "\n";
      printGrid(cur);
    }

    // Try every possible single-position change
    for (int pos = 0; pos < g_numPos; pos++) {
      int curType = stateGet(cur, pos);
      for (int newType = 0; newType <= g_maxType[pos]; newType++) {
        if (newType == curType)
          continue;
        if (!canPlace(cur, pos, newType))
          continue;

        uint64_t next = stateSet(cur, pos, newType);
        int nextScore = curScore + (newType - curType);
        if (visited.insert({next, {cur, (int8_t)pos, (int8_t)newType}})
                .second) {
          q.push_back({next, nextScore});
        }
      }
    }

    if (explored % 2000000 == 0)
      cout << "Explored: " << explored << " | Visited: " << visited.size()
           << " | Best: " << bestScore << "\r" << flush;
  }

  cout << "\nBFS complete.\n";
  cout << "States explored : " << explored << "\n";
  cout << "Unique states   : " << visited.size() << "\n";
  cout << "Best score      : " << bestScore << "\n\n";
  cout << "Optimal layout:\n";
  printGrid(bestState);

  printBuildPath(bestState, visited);

  return 0;
}
