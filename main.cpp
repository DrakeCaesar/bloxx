#include <algorithm>
#include <iostream>
#include <map>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

// ANSI color codes
const string COLORS[] = {
    "\033[44m", // Blue background   - type 0
    "\033[41m", // Red background    - type 1
    "\033[42m", // Green background  - type 2
    "\033[43m", // Yellow background - type 3
    "\033[47m"  // White background  - type 4
};
const string COLOR_NAMES[] = {"Blue", "Red", "Green", "Yellow", "White"};
const string RESET = "\033[0m";
const string BLACK_TEXT = "\033[30m"; // Black text for number

// ══════════════════════════════════════════════════════════════════════════
// Grid Configuration
// ══════════════════════════════════════════════════════════════════════════

struct GridConfig {
  int rows;
  int cols;
  vector<pair<int, int>> positions;    // list of (row, col) coords
  map<pair<int, int>, int> coordToIdx; // coord -> position index
  vector<vector<int>> neighbors; // neighbors[idx] = list of neighbor indices
  vector<int> maxType;           // maxType[idx] = max tower type achievable

  int numPositions() const { return (int)positions.size(); }

  void computeNeighbors() {
    neighbors.clear();
    neighbors.resize(numPositions());

    for (int i = 0; i < numPositions(); i++) {
      int r = positions[i].first;
      int c = positions[i].second;
      // Check 4 orthogonal directions
      vector<pair<int, int>> dirs = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
      for (int d = 0; d < 4; d++) {
        int dr = dirs[d].first;
        int dc = dirs[d].second;
        pair<int, int> neighbor = {r + dr, c + dc};
        if (coordToIdx.count(neighbor)) {
          neighbors[i].push_back(coordToIdx[neighbor]);
        }
      }
    }
  }

  void computeMaxTypes() {
    maxType.clear();
    maxType.resize(numPositions());
    for (int i = 0; i < numPositions(); i++) {
      maxType[i] = min(4, (int)neighbors[i].size());
    }
  }

  void build() {
    coordToIdx.clear();
    for (int i = 0; i < numPositions(); i++) {
      coordToIdx[positions[i]] = i;
    }
    computeNeighbors();
    computeMaxTypes();
  }
};

// Default: 4x4 grid without corners (Tower Blox Deluxe first level)
GridConfig createDefaultGrid() {
  GridConfig g;
  g.rows = 4;
  g.cols = 4;
  g.positions = {
      {0, 1}, {0, 2},                 // row 0: positions 0, 1
      {1, 0}, {1, 1}, {1, 2}, {1, 3}, // row 1: positions 2, 3, 4, 5
      {2, 0}, {2, 1}, {2, 2}, {2, 3}, // row 2: positions 6, 7, 8, 9
      {3, 1}, {3, 2}                  // row 3: positions 10, 11
  };
  g.build();
  return g;
}

GridConfig GRID = createDefaultGrid();

// Create grid from string pattern
// '.' = valid position, ' ' or '#' = blank/invalid
GridConfig createGridFromPattern(const vector<string> &pattern) {
  GridConfig g;
  g.rows = (int)pattern.size();
  g.cols = 0;
  for (const auto &row : pattern) {
    g.cols = max(g.cols, (int)row.size());
  }

  for (int r = 0; r < g.rows; r++) {
    for (int c = 0; c < (int)pattern[r].size(); c++) {
      if (pattern[r][c] == '.') {
        g.positions.push_back({r, c});
      }
    }
  }

  g.build();
  return g;
}

// Example patterns
GridConfig createSmallGrid() { return createGridFromPattern({"..", ".."}); }

GridConfig createCrossGrid() {
  return createGridFromPattern({" . ", "...", " . "});
}

GridConfig createLargeGrid() {
  return createGridFromPattern({".....", ".....", ".....", ".....", "....."});
}

GridConfig createCustomGrid() {
  return createGridFromPattern({" ... ", ".....", ".....", ".....", " ... "});
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

// ── Grid display ───────────────────────────────────────────────────────────

void printGrid(uint64_t s) {
  auto cell = [&](int r, int c) -> string {
    pair<int, int> coord = {r, c};
    if (!GRID.coordToIdx.count(coord))
      return " ";
    int idx = GRID.coordToIdx[coord];
    int t = stateGet(s, idx);
    return COLORS[t] + BLACK_TEXT + to_string(t) + RESET;
  };

  for (int r = 0; r < GRID.rows; r++) {
    cout << "  ";
    for (int c = 0; c < GRID.cols; c++) {
      cout << cell(r, c);
      if (c < GRID.cols - 1)
        cout << " ";
    }
    cout << "\n";
  }
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
         << BLACK_TEXT << COLOR_NAMES[oldType] << RESET << " (" << oldType
         << ") -> " << COLORS[newType] << BLACK_TEXT << COLOR_NAMES[newType]
         << RESET << " (" << newType << ")\n";
    printGrid(state);
  }
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
