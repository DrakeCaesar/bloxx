#pragma once

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
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

// Create grid from string pattern
// '.' = valid position, ' ' or '#' = blank/invalid
inline GridConfig createGridFromPattern(const vector<string> &pattern) {
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

// Game grids (Tower Blox Deluxe levels)
inline GridConfig createLevel1Grid() {
  return createGridFromPattern({" .. ", "....", "....", " .. "});
}

inline GridConfig createLevel2Grid() {
  return createGridFromPattern({"....", "....", "....", "...."});
}

inline GridConfig createLevel3Grid() {
  return createGridFromPattern({".....", ".....", "  .  ", ".....", "....."});
}

inline GridConfig createLevel4Grid() {
  return createGridFromPattern(
      {"......", "..  ..", " .... ", " .... ", "..  ..", "......"});
};

inline GridConfig createLevel5Grid() {
  return createGridFromPattern({"    ...", " ......", " ......", " .. .. ",
                                "...... ", "...... ", "...    "});
};

const int LEVEL = 2;
// 0 = BFS (queue), 1 = DFS (stack), 2 = IDDFS (iterative deepening DFS)
const int SEARCH_MODE = 1;

inline GridConfig createGrid() {
  switch (LEVEL) {
  case 2:
    return createLevel2Grid();
  case 3:
    return createLevel3Grid();
  case 4:
    return createLevel4Grid();
  case 5:
    return createLevel5Grid();
  default:
    return createLevel1Grid();
  }
}

inline GridConfig &getGrid() {
  static GridConfig GRID = createGrid();
  return GRID;
}

#define GRID getGrid()

// ── Grid display ───────────────────────────────────────────────────────────

// Forward declarations for state functions (defined in main.cpp)
extern int stateGet(uint64_t s, int pos);
extern uint64_t stateSet(uint64_t s, int pos, int val);

inline void printGrid(uint64_t s) {
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
// printBuildPath is defined in main.cpp (needs canonical(), which lives there)
