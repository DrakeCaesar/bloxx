#include <algorithm>
#include <iostream>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

// ANSI color codes
const string COLORS[] = {
    "\033[34m", // Blue   - type 0
    "\033[31m", // Red    - type 1
    "\033[32m", // Green  - type 2
    "\033[33m", // Yellow - type 3
    "\033[37m"  // White  - type 4
};
const string COLOR_NAMES[] = {"Blue", "Red", "Green", "Yellow", "White"};
const string RESET = "\033[0m";

// Grid layout (orthogonal adjacency only):
//   X  0  1  X
//   2  3  4  5
//   6  7  8  9
//   X  10 11 X

const vector<int> NEIGHBORS[12] = {
    {1, 3},        // 0
    {0, 4},        // 1
    {3, 6},        // 2
    {0, 2, 4, 7},  // 3
    {1, 3, 5, 8},  // 4
    {4, 9},        // 5
    {2, 7},        // 6
    {3, 6, 8, 10}, // 7
    {4, 7, 9, 11}, // 8
    {5, 8},        // 9
    {7, 11},       // 10
    {8, 10}        // 11
};

// Max achievable type = min(neighbor count, 4)
const int MAX_TYPE[12] = {
    2, // 0:  2 neighbors
    2, // 1:  2 neighbors
    2, // 2:  2 neighbors
    4, // 3:  4 neighbors
    4, // 4:  4 neighbors
    2, // 5:  2 neighbors
    2, // 6:  2 neighbors
    4, // 7:  4 neighbors
    4, // 8:  4 neighbors
    2, // 9:  2 neighbors
    2, // 10: 2 neighbors
    2  // 11: 2 neighbors
};

// ── State packing ──────────────────────────────────────────────────────────
// Each position stores a value 0-4 (fits in 3 bits).
// 12 positions × 3 bits = 36 bits → fits in uint64_t.

inline uint64_t stateSet(uint64_t s, int pos, int val) {
  s &= ~(uint64_t(7) << (pos * 3));
  s |= uint64_t(val) << (pos * 3);
  return s;
}
inline int stateGet(uint64_t s, int pos) { return (s >> (pos * 3)) & 7; }
inline int stateScore(uint64_t s) {
  int total = 0;
  for (int i = 0; i < 12; i++)
    total += stateGet(s, i);
  return total;
}

// ── Placement rule ─────────────────────────────────────────────────────────
// Tower of type T at pos requires one neighbor of each type 0..T-1.

bool canPlace(uint64_t state, int pos, int type) {
  if (type == 0)
    return true;
  bool has[5] = {};
  for (int n : NEIGHBORS[pos])
    has[stateGet(state, n)] = true;
  for (int i = 0; i < type; i++)
    if (!has[i])
      return false;
  return true;
}

// ── Grid display ───────────────────────────────────────────────────────────

void printGrid(uint64_t s) {
  auto cell = [&](int pos) {
    int t = stateGet(s, pos);
    return COLORS[t] + to_string(t) + RESET;
  };
  cout << "  X  " << cell(0) << "  " << cell(1) << "  X\n";
  cout << "  ";
  for (int i = 2; i <= 5; i++)
    cout << cell(i) << "  ";
  cout << "\n  ";
  for (int i = 6; i <= 9; i++)
    cout << cell(i) << "  ";
  cout << "\n  X  " << cell(10) << "  " << cell(11) << "  X\n";
}

// ── Path reconstruction ────────────────────────────────────────────────────

struct Transition {
  uint64_t parent;
  int8_t pos;  // position changed (-1 = start sentinel)
  int8_t type; // new type applied
};

void printBuildPath(uint64_t bestState,
                    const unordered_map<uint64_t, Transition> &visited) {
  // Walk back from best state to start
  vector<pair<int, int>> steps; // (pos, newType)
  uint64_t cur = bestState;
  while (visited.at(cur).pos != -1) {
    const auto &tr = visited.at(cur);
    steps.push_back({tr.pos, tr.type});
    cur = tr.parent;
  }
  reverse(steps.begin(), steps.end());

  cout << "\n--- Build path (" << steps.size() << " steps) ---\n\n";

  uint64_t state = 0; // start: all type 0
  cout << "Start:\n";
  printGrid(state);

  for (int i = 0; i < (int)steps.size(); i++) {
    int pos = steps[i].first;
    int newType = steps[i].second;
    int oldType = stateGet(state, pos);
    state = stateSet(state, pos, newType);

    cout << "\nStep " << (i + 1) << ": pos " << pos << "  " << COLORS[oldType]
         << COLOR_NAMES[oldType] << RESET << " (" << oldType << ") -> "
         << COLORS[newType] << COLOR_NAMES[newType] << RESET << " (" << newType
         << ")\n";
    printGrid(state);
  }
}

// ── Main ───────────────────────────────────────────────────────────────────

int main() {
  cout << "Tower Blox Deluxe Solver\n";
  cout << "========================\n";
  cout << "Tower types:\n";
  for (int i = 0; i <= 4; i++)
    cout << "  " << COLORS[i] << i << RESET << " (" << COLOR_NAMES[i] << ")\n";
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
    for (int pos = 0; pos < 12; pos++) {
      int curType = stateGet(cur, pos);
      for (int newType = 0; newType <= MAX_TYPE[pos]; newType++) {
        if (newType == curType)
          continue;
        if (!canPlace(cur, pos, newType))
          continue;

        uint64_t next = stateSet(cur, pos, newType);
        if (!visited.count(next)) {
          visited[next] = {cur, (int8_t)pos, (int8_t)newType};
          q.push(next);
        }
      }
    }

    if (explored % 200000 == 0)
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
