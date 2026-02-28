#pragma once

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

// ════════════════════════════════════════════════════════════════════════════
// State Type Selection (based on grid size)
// ════════════════════════════════════════════════════════════════════════════
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

// ════════════════════════════════════════════════════════════════════════════
// Configuration Constants
// ════════════════════════════════════════════════════════════════════════════

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
    // 3: 5x5 square (25 positions)
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

// ANSI color codes and display configuration
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

