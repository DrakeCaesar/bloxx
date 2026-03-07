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
// Each position needs 3 bits (types 0-5, 0=empty). 128-bit → 42 positions max.

// Windows (MSVC or Clang): 128-bit state in __m128i (SSE2/AVX2), passed in XMM
// for stable ABI.
#include <immintrin.h>
struct uint128_simd {
  union {
    __m128i v;
    struct {
      uint64_t lo;
      uint64_t hi;
    } u;
  };
  uint128_simd() : v(_mm_setzero_si128()) {}
  uint128_simd(uint64_t l, uint64_t h = 0) {
    u.lo = l;
    u.hi = h;
  }
  uint128_simd(int x) {
    u.lo = static_cast<uint64_t>(x);
    u.hi = 0;
  }
  uint128_simd operator<<(int n) const {
    if (n <= 0)
      return *this;
    if (n >= 64)
      return uint128_simd(0, u.lo << (n - 64));
    return uint128_simd(u.lo << n, (u.hi << n) | (u.lo >> (64 - n)));
  }
  uint128_simd operator>>(int n) const {
    if (n <= 0)
      return *this;
    if (n >= 64)
      return uint128_simd(u.hi >> (n - 64), 0);
    return uint128_simd((u.lo >> n) | (u.hi << (64 - n)), u.hi >> n);
  }
  uint128_simd operator&(uint128_simd o) const {
    uint128_simd r;
    r.v = _mm_and_si128(v, o.v);
    return r;
  }
  uint128_simd operator&(uint64_t m) const { return uint128_simd(u.lo & m, 0); }
  uint128_simd operator&(int m) const {
    return uint128_simd(u.lo & static_cast<uint64_t>(m), 0);
  }
  uint128_simd operator|(uint128_simd o) const {
    uint128_simd r;
    r.v = _mm_or_si128(v, o.v);
    return r;
  }
  uint128_simd operator~() const {
    uint128_simd r;
    r.v = _mm_xor_si128(v, _mm_set1_epi64x(-1));
    return r;
  }
  uint128_simd &operator&=(uint128_simd o) {
    v = _mm_and_si128(v, o.v);
    return *this;
  }
  uint128_simd &operator|=(uint128_simd o) {
    v = _mm_or_si128(v, o.v);
    return *this;
  }
  bool operator==(uint128_simd o) const {
    return u.lo == o.u.lo && u.hi == o.u.hi;
  }
  bool operator==(int o) const {
    return u.hi == 0 && u.lo == static_cast<uint64_t>(o);
  }
  bool operator!=(uint128_simd o) const {
    return u.lo != o.u.lo || u.hi != o.u.hi;
  }
  bool operator!=(int o) const {
    return u.hi != 0 || u.lo != static_cast<uint64_t>(o);
  }
  bool operator<(uint128_simd o) const {
    return u.hi != o.u.hi ? u.hi < o.u.hi : u.lo < o.u.lo;
  }
  operator uint64_t() const { return u.lo; }
  operator int() const { return static_cast<int>(u.lo); }
};
namespace std {
template <> struct hash<uint128_simd> {
  size_t operator()(const uint128_simd &x) const noexcept {
    return std::hash<uint64_t>{}(x.u.lo) ^ (std::hash<uint64_t>{}(x.u.hi) << 1);
  }
};
} // namespace std

using StateType64 = uint64_t;
using StateType128 = uint128_simd;

// ════════════════════════════════════════════════════════════════════════════
// Configuration Constants
// ════════════════════════════════════════════════════════════════════════════

#define LEVEL 2 // Grid level (1-5), must be #define for preprocessor
// 0 = BFS, 1 = Sliding Window
[[maybe_unused]] const int SEARCH_MODE = 1;
// Maximum tower type unlocked (1=Blue, 2=Red, 3=Green, 4=Yellow, 5=White)
// Set to 5 for all towers unlocked, 3 for only Blue/Red/Green, etc.
const int MAX_UNLOCKED_TYPE = 5;

// Baseline generation: generate baselines for unlocked types in this range
const int UNLOCKED_TYPE_START = 2; // First type to generate baseline for
const int UNLOCKED_TYPE_END = 5;   // Last type to generate baseline for

// Multi-level configuration
const int LEVEL_START = 1; // First level to process
const int LEVEL_END = 1;   // Last level to process (inclusive)

// Current level being processed (for multi-level runs)
inline int &getCurrentLevel() {
  static int currentLevel = LEVEL;
  return currentLevel;
}

// Reverse BFS (path reconstruction): number of threads per layer
const int REVERSE_BFS_NUM_THREADS = 8;

// When true, path-finding bidirectional BFS starts from empty (all 0) instead
// of all-blue (all 1). Can find more optimal paths (e.g. fewer placements).
const bool PATH_FROM_EMPTY = true;

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
// Index 0 = empty, 1-5 = tower types
const string COLORS[] = {
    "\033[100m", // Grey background   - type 0 (empty)
    "\033[44m",  // Blue background   - type 1
    "\033[41m",  // Red background    - type 2
    "\033[42m",  // Green background  - type 3
    "\033[43m",  // Yellow background - type 4
    "\033[47m"   // White background  - type 5
};
const string COLOR_NAMES[] = {"Empty", "Blue",   "Red",
                              "Green", "Yellow", "White"};
const string RESET = "\033[0m";
const string BLACK_TEXT = "\033[30m"; // Black text for number
const string EMPTY_BG =
    "\033[100m"; // Grey background - empty/unplaced position
const string EMPTY_FG =
    "\033[90m"; // Grey foreground - empty position (grey bg + grey fg)

// ANSI foreground colors for half-block separators (3x codes)
// Index 0 = empty, 1-5 = tower types
const string FG_COLORS[] = {
    "\033[90m", // Grey foreground   - type 0 (empty)
    "\033[34m", // Blue foreground   - type 1
    "\033[31m", // Red foreground    - type 2
    "\033[32m", // Green foreground  - type 3
    "\033[33m", // Yellow foreground - type 4
    "\033[37m"  // White foreground  - type 5
};
