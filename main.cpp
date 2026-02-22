#include <iostream>
#include <queue>
#include <unordered_map>
#include <vector>

#include "grid.h"

using namespace std;

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
  for (int i = 0; i < GRID.numPositions(); i++)
    total += stateGet(s, i);
  return total;
}

// ── Placement rule ─────────────────────────────────────────────────────────
// Tower of type T at pos requires one neighbor of each type 0..T-1.

bool canPlace(uint64_t state, int pos, int type) {
  if (type == 0)
    return true;
  bool has[5] = {};
  for (int n : GRID.neighbors[pos])
    has[stateGet(state, n)] = true;
  for (int i = 0; i < type; i++)
    if (!has[i])
      return false;
  return true;
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

  // BFS: every state is visited exactly once
  unordered_map<uint64_t, Transition> visited;
  visited.reserve(1 << 23); // pre-allocate ~8M buckets
  queue<uint64_t> q;

  uint64_t startState = 0;
  visited[startState] = {0, -1, -1}; // sentinel
  q.push(startState);

  uint64_t bestState = startState;
  int bestScore = 0;
  long long explored = 0;

  while (!q.empty()) {
    uint64_t cur = q.front();
    q.pop();
    explored++;

    int score = stateScore(cur);
    if (score > bestScore) {
      bestScore = score;
      bestState = cur;
      cout << "New best score: " << bestScore << "\n";
      printGrid(cur);
    }

    // Try every possible single-position change
    for (int pos = 0; pos < GRID.numPositions(); pos++) {
      int curType = stateGet(cur, pos);
      for (int newType = 0; newType <= GRID.maxType[pos]; newType++) {
        if (newType == curType)
          continue;
        if (!canPlace(cur, pos, newType))
          continue;

        uint64_t next = stateSet(cur, pos, newType);
        if (visited.insert({next, {cur, (int8_t)pos, (int8_t)newType}})
                .second) {
          q.push(next);
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
