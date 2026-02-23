#pragma once

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

// ══════════════════════════════════════════════════════════════════════════
// State Type Selection (based on grid size)
// ══════════════════════════════════════════════════════════════════════════
// Each position needs 3 bits (types 0-4)
// uint64_t: 64 bits / 3 = 21 positions max
// Custom 128-bit using SIMD: 128 bits / 3 = 42 positions max

// Use native __uint128_t (GCC/Clang compiler built-in)
using uint128_simd = __uint128_t;

using StateType64 = uint64_t;
using StateType128 = uint128_simd;

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

  void computeMaxTypes(int maxUnlockedType = 4) {
    maxType.clear();
    maxType.resize(numPositions());
    for (int i = 0; i < numPositions(); i++) {
      // Max type is limited by: neighbors available AND towers unlocked
      maxType[i] = min(maxUnlockedType, min(4, (int)neighbors[i].size()));
    }
  }

  void build(int maxUnlockedType = 4) {
    coordToIdx.clear();
    for (int i = 0; i < numPositions(); i++) {
      coordToIdx[positions[i]] = i;
    }
    computeNeighbors();
    computeMaxTypes(maxUnlockedType);
  }
};

// ══════════════════════════════════════════════════════════════════════════
// Configuration Constants
// ══════════════════════════════════════════════════════════════════════════

#define LEVEL 2  // Grid level (1-5), must be #define for preprocessor
// 0 = BFS (queue), 1 = DFS (stack), 2 = IDDFS (iterative deepening DFS), 3 = Sliding Window
[[maybe_unused]] const int SEARCH_MODE = 3;
// Maximum tower type unlocked (0=Blue, 1=Red, 2=Green, 3=Yellow, 4=White)
// Set to 4 for all towers unlocked, 2 for only Blue/Red/Green, etc.
const int MAX_UNLOCKED_TYPE = 4;
// Sliding window parameters
const int MAX_DEPTH = 10;  // Maximum BFS depth per window

// Kernel patterns: '.' = include in window, ' ' = exclude
// Pattern is centered on the window position
// Multiple kernels can be defined; select one with KERNEL_INDEX
const vector<vector<string>> KERNEL_PATTERNS = {
    // 0: Cross shape (5 positions)
    {" . ", "...", " . "},
    // 1: 3x3 square (9 positions)
    {"...", "...", "..."},
    // 2: Diamond / plus shape (13 positions)
    {"  .  ", " ... ", ".....", " ... ", "  .  "},
    // 3: 5x5 square (25 positions)
    {".....", ".....", ".....", ".....", "....."},
    // 4: Plus shape (9 positions)
    {"  .  ", "  .  ", ".....", "  .  ", "  .  "},
    // 5: L-shape (5 positions)
    {"...", ". ", ". "},
};
const int KERNEL_INDEX = 2;  // Select which kernel to use
const vector<string> &WINDOW_KERNEL_PATTERN = KERNEL_PATTERNS[KERNEL_INDEX];

// Select state type based on level - Level 1-3 use 64-bit, Level 4-5 use
// 128-bit
#if LEVEL <= 3
using StateType = StateType64;
#else
using StateType = StateType128;
#endif

// ══════════════════════════════════════════════════════════════════════════
// Grid Creation
// ══════════════════════════════════════════════════════════════════════════

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

  g.build(MAX_UNLOCKED_TYPE);
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
extern int stateGet(StateType s, int pos);
extern StateType stateSet(StateType s, int pos, int val);

inline void printGrid(StateType s) {
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

// Transition struct optimized for board size
// For grids ≤52 positions: pack pos (6 bits) + type (3 bits) + sentinel into
// uint16_t This reduces Transition size and improves cache utilization
struct Transition {
  StateType parent;
  uint16_t packed; // bits 0-5: pos (0-63, 63=sentinel), bits 6-8: type (0-4)

  Transition(StateType p = 0, int8_t pos = -1, int8_t type = 0)
      : parent(p), packed((pos == -1 ? 63 : pos) | (type << 6)) {}

  int8_t getPos() const { return (packed & 0x3F) == 63 ? -1 : (packed & 0x3F); }
  int8_t getType() const { return (packed >> 6) & 0x7; }
};

// printBuildPath is defined in main.cpp (needs canonical(), which lives there)
