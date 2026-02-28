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

// Hash specialization for __uint128_t (needed for robin_hood and std maps)
// Only define if not already provided by standard library (GCC 13+ has it)
#if !defined(__GLIBCXX__) || (__GLIBCXX__ < 20230528)
namespace std {
template <> struct hash<__uint128_t> {
  size_t operator()(__uint128_t x) const noexcept {
    // Simple hash: XOR the upper and lower 64 bits
    uint64_t low = static_cast<uint64_t>(x);
    uint64_t high = static_cast<uint64_t>(x >> 64);
    return std::hash<uint64_t>{}(low) ^ (std::hash<uint64_t>{}(high) << 1);
  }
};
} // namespace std
#endif

// ANSI color codes

// ══════════════════════════════════════════════════════════════════════════
// Configuration Constants
// ══════════════════════════════════════════════════════════════════════════

#define LEVEL 2 // Grid level (1-5), must be #define for preprocessor
// 0 = BFS (queue), 1 = DFS (stack), 2 = IDDFS (iterative deepening DFS), 3 =
// Sliding Window
[[maybe_unused]] const int SEARCH_MODE = 3;
// Maximum tower type unlocked (0=Blue, 1=Red, 2=Green, 3=Yellow, 4=White)
// Set to 4 for all towers unlocked, 2 for only Blue/Red/Green, etc.
const int MAX_UNLOCKED_TYPE = 4;

// Baseline generation: generate baselines for unlocked types in this range
const int UNLOCKED_TYPE_START = 1; // First type to generate baseline for
const int UNLOCKED_TYPE_END = 4;   // Last type to generate baseline for

// Multi-level configuration
const int LEVEL_START = 1; // First level to process
const int LEVEL_END = 5;   // Last level to process (inclusive)

// Current level being processed (for multi-level runs)
inline int &getCurrentLevel() {
  static int currentLevel = LEVEL;
  return currentLevel;
}
// Sliding window parameters
const int MAX_DEPTH = 10; // Maximum BFS depth per window
const int SCAFFOLD_DEPTH =
    10; // Maximum BFS depth for scaffold search during path reconstruction

const int KERNEL_INDEX = 2; // Select which kernel to use

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
    // 3: 5x5 square (25 positiogitns)
    {".....", ".....", ".....", ".....", "....."},
    // 4: Plus shape (9 positions)
    {"  .  ", "  .  ", ".....", "  .  ", "  .  "},
    // 5: L-shape (5 positions)
    {"...", ". ", ". "},
};
inline const vector<string> &getWindowKernelPattern() {
  static const vector<string> pattern = KERNEL_PATTERNS[KERNEL_INDEX];
  return pattern;
}
#define WINDOW_KERNEL_PATTERN getWindowKernelPattern()

// Select state type based on level - Level 1-3 use 64-bit, Level 4-5 use
// 128-bit
// For multi-level processing, we use 128-bit for all levels to avoid
// recompilation Comment out the #if when running single level, use the forced
// line for multi-level
/*
#if LEVEL <= 3
using StateType = StateType64;
#else
using StateType = StateType128;
#endif
*/

// Force 128-bit for multi-level support
using StateType = StateType128;

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
const string EMPTY_BG =
    "\033[100m"; // Grey background - empty/unplaced position
const string EMPTY_FG =
    "\033[90m"; // Grey foreground - empty position (grey bg + grey fg)

// ANSI foreground colors for half-block separators (3x codes)
const string FG_COLORS[] = {
    "\033[34m", // Blue foreground   - type 0
    "\033[31m", // Red foreground    - type 1
    "\033[32m", // Green foreground  - type 2
    "\033[33m", // Yellow foreground - type 3
    "\033[37m"  // White foreground  - type 4
};

inline const string &bgColor(int t, bool placed) {
  return placed ? COLORS[t] : EMPTY_BG;
}
inline const string &fgColor(int t, bool placed) {
  return placed ? FG_COLORS[t] : EMPTY_FG;
}

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

inline GridConfig createGrid(int level) {
  switch (level) {
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

inline void reinitializeGrid(int level) {
  getGrid() = createGrid(level);
  getCurrentLevel() = level;
}

#define GRID getGrid()

// ── Grid display ───────────────────────────────────────────────────────────

// Forward declarations for state functions (defined in main.cpp)
extern int stateGet(StateType s, int pos);
extern StateType stateSet(StateType s, int pos, int val);

inline void printGrid(StateType s) {
  for (int r = 0; r < GRID.rows; r++) {
    pair<int, int> firstCoord = {r, 0};
    if (GRID.coordToIdx.count(firstCoord)) {
      int firstT = stateGet(s, GRID.coordToIdx[firstCoord]);
      cout << FG_COLORS[firstT] << "▐" << RESET;
    } else {
      cout << " ";
    }

    for (int c = 0; c < GRID.cols; c++) {
      pair<int, int> coord = {r, c};
      if (GRID.coordToIdx.count(coord)) {
        int t = stateGet(s, GRID.coordToIdx[coord]);
        cout << COLORS[t] << BLACK_TEXT << t << RESET;
      } else {
        cout << " ";
      }

      if (c < GRID.cols - 1) {
        pair<int, int> nextCoord = {r, c + 1};
        bool currExists = GRID.coordToIdx.count(coord);
        bool nextExists = GRID.coordToIdx.count(nextCoord);

        if (currExists && nextExists) {
          int t = stateGet(s, GRID.coordToIdx[coord]);
          int nextT = stateGet(s, GRID.coordToIdx[nextCoord]);
          cout << FG_COLORS[nextT] << COLORS[t] << "▐" << RESET;
        } else if (currExists && !nextExists) {
          int t = stateGet(s, GRID.coordToIdx[coord]);
          cout << FG_COLORS[t] << "▌" << RESET;
        } else if (!currExists && nextExists) {
          int nextT = stateGet(s, GRID.coordToIdx[nextCoord]);
          cout << FG_COLORS[nextT] << "▐" << RESET;
        } else {
          cout << " ";
        }
      }
    }

    pair<int, int> lastCoord = {r, GRID.cols - 1};
    if (GRID.coordToIdx.count(lastCoord)) {
      int lastT = stateGet(s, GRID.coordToIdx[lastCoord]);
      cout << FG_COLORS[lastT] << "▌" << RESET << "\n";
    } else {
      cout << " \n";
    }
  }
}

// printGrid with placed mask: bit i=1 means position i has been explicitly
// built. Unplaced positions (bit=0, type=0) are rendered as dark empty cells.
inline void printGrid(StateType s, uint64_t placedMask) {
  for (int r = 0; r < GRID.rows; r++) {
    pair<int, int> firstCoord = {r, 0};
    if (GRID.coordToIdx.count(firstCoord)) {
      int idx = GRID.coordToIdx[firstCoord];
      int firstT = stateGet(s, idx);
      bool p = (placedMask >> idx) & 1;
      cout << fgColor(firstT, p) << "▐" << RESET;
    } else {
      cout << " ";
    }

    for (int c = 0; c < GRID.cols; c++) {
      pair<int, int> coord = {r, c};
      if (GRID.coordToIdx.count(coord)) {
        int idx = GRID.coordToIdx[coord];
        int t = stateGet(s, idx);
        bool p = (placedMask >> idx) & 1;
        if (p) {
          cout << COLORS[t] << BLACK_TEXT << t << RESET;
        } else {
          cout << EMPTY_BG << EMPTY_FG << "·" << RESET;
        }
      } else {
        cout << " ";
      }

      if (c < GRID.cols - 1) {
        pair<int, int> nextCoord = {r, c + 1};
        bool currExists = GRID.coordToIdx.count(coord);
        bool nextExists = GRID.coordToIdx.count(nextCoord);

        if (currExists && nextExists) {
          int idx = GRID.coordToIdx[coord];
          int nextIdx = GRID.coordToIdx[nextCoord];
          int t = stateGet(s, idx);
          int nextT = stateGet(s, nextIdx);
          bool p = (placedMask >> idx) & 1;
          bool np = (placedMask >> nextIdx) & 1;
          cout << fgColor(nextT, np) << bgColor(t, p) << "▐" << RESET;
        } else if (currExists && !nextExists) {
          int idx = GRID.coordToIdx[coord];
          int t = stateGet(s, idx);
          bool p = (placedMask >> idx) & 1;
          cout << fgColor(t, p) << "▌" << RESET;
        } else if (!currExists && nextExists) {
          int nextIdx = GRID.coordToIdx[nextCoord];
          int nextT = stateGet(s, nextIdx);
          bool np = (placedMask >> nextIdx) & 1;
          cout << fgColor(nextT, np) << "▐" << RESET;
        } else {
          cout << " ";
        }
      }
    }

    pair<int, int> lastCoord = {r, GRID.cols - 1};
    if (GRID.coordToIdx.count(lastCoord)) {
      int idx = GRID.coordToIdx[lastCoord];
      int lastT = stateGet(s, idx);
      bool p = (placedMask >> idx) & 1;
      cout << fgColor(lastT, p) << "▌" << RESET << "\n";
    } else {
      cout << " \n";
    }
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
