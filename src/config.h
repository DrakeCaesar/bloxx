#pragma once

#include <cstddef>
#include <cstdint>
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
  uint128_simd operator&(int m) const { return uint128_simd(u.lo & static_cast<uint64_t>(m), 0); }
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
  bool operator==(uint128_simd o) const { return u.lo == o.u.lo && u.hi == o.u.hi; }
  bool operator==(int o) const { return u.hi == 0 && u.lo == static_cast<uint64_t>(o); }
  bool operator!=(uint128_simd o) const { return u.lo != o.u.lo || u.hi != o.u.hi; }
  bool operator!=(int o) const { return u.hi != 0 || u.lo != static_cast<uint64_t>(o); }
  bool operator<(uint128_simd o) const { return u.hi != o.u.hi ? u.hi < o.u.hi : u.lo < o.u.lo; }
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

// 0 = BFS, 1 = Sliding Window, 2 = DFS (depth-limited, visited set + path
// stack)
[[maybe_unused]] const int SEARCH_MODE = 1;
// DFS mode only: max depth for depth-limited DFS
const int DFS_MAX_DEPTH = 17;

// Baseline generation: generate baselines for unlocked types in this range
const int UNLOCKED_TYPE_START = 5; // First type to generate baseline for
const int UNLOCKED_TYPE_END = 5;   // Last type to generate baseline for

// Multi-level configuration
const int LEVEL_START = 2; // First level to process
const int LEVEL_END = 2;   // Last level to process (inclusive)

// Current level being processed (for multi-level runs)
inline int &getCurrentLevel() {
  static int currentLevel = LEVEL_START;
  return currentLevel;
}

// Reverse BFS (path reconstruction): number of threads per layer
const int REVERSE_BFS_NUM_THREADS = 4;
// Throttle for status-line updates (expand/prune/merge progress); milliseconds
const int REVERSE_BFS_STATUS_UPDATE_MS = 1000;

// Reverse BFS pruning: only run when parent map size exceeds this
const size_t REVERSE_BFS_PRUNE_THRESHOLD = 1'000'000;
// Reverse side: keep states with tower count <= minInQueue + this delta
const int REVERSE_BFS_PRUNE_LOW_DELTA = 3;
// Forward side: keep states with tower count >= maxInQueue - this delta
const int REVERSE_BFS_PRUNE_HIGH_DELTA = 3;

// When true, run forward BFS in parallel with reverse and meet in the middle
// (finds shorter path, uses more memory). When false, only reverse BFS runs;
// path is built by walking from final back to start (less memory, path may be
// longer).
const bool RUN_PARALLEL_FORWARD_BFS = false;

// When true, path-finding bidirectional BFS starts from empty (all 0) instead
// of all-blue (all 1). Can find more optimal paths (e.g. fewer placements).
const bool PATH_FROM_EMPTY = true;

// ── Transition tables ────────────────────────────────────────────────────────
// TRANSITIONS[fromType][toType] = true if (fromType → toType) is allowed.
// Edit these to prune transitions that are never optimal.
// The canPlace() neighbor check is still applied on top of these.
//                empty  blue   red    green  yellow white
// Normal forward BFS (and bidirectional forward): no downgrading from white.
constexpr bool FORWARD_TRANSITIONS[6][6] = {
    /* empty  */ {0, 1, 1, 1, 1, 1},
    /* blue   */ {0, 0, 1, 1, 1, 1},
    /* red    */ {0, 1, 0, 1, 1, 1},
    /* green  */ {0, 1, 1, 0, 1, 1},
    /* yellow */ {0, 1, 1, 1, 0, 1},
    /* white  */ {0, 0, 0, 0, 0, 0},
};

// Sliding window only: allow white→1,2,3,4 so depth-limited BFS can tear down
// and rebuild for better local optima.
constexpr bool FORWARD_TRANSITIONS_SLIDING_WINDOW[6][6] = {
    /* empty  */ {0, 1, 1, 1, 1, 1},
    /* blue   */ {0, 0, 1, 1, 1, 1},
    /* red    */ {0, 1, 0, 1, 1, 1},
    /* green  */ {0, 1, 1, 0, 1, 1},
    /* yellow */ {0, 1, 1, 1, 0, 1},
    /* white  */ {0, 1, 1, 1, 1, 0},
};

constexpr bool BACKWARD_TRANSITIONS[6][6] = {
    /* empty  */ {0, 0, 0, 0, 0, 0},
    /* blue   */ {0, 0, 1, 1, 1, 0},
    /* red    */ {0, 1, 0, 1, 1, 0},
    /* green  */ {0, 1, 1, 0, 1, 0},
    /* yellow */ {0, 1, 1, 1, 0, 0},
    /* white  */ {0, 1, 1, 1, 1, 0},
};

constexpr bool BACKWARD_TRANSITIONS_EMPTY[6][6] = {
    /* empty  */ {0, 0, 0, 0, 0, 0},
    /* blue   */ {1, 0, 1, 1, 1, 0},
    /* red    */ {1, 1, 0, 1, 1, 0},
    /* green  */ {1, 1, 1, 0, 1, 0},
    /* yellow */ {1, 1, 1, 1, 0, 0},
    /* white  */ {1, 1, 1, 1, 1, 0},
};

// Sliding window parameters
const int MAX_DEPTH = 10; // Maximum BFS depth per window

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
// 128-bit. For multi-level we use 128-bit for all levels. For single-level
// you could switch at runtime using getCurrentLevel() <= 3.
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
const string COLOR_NAMES[] = {"Empty", "Blue", "Red", "Green", "Yellow", "White"};
const string RESET = "\033[0m";
const string BLACK_TEXT = "\033[30m"; // Black text for number
const string EMPTY_BG = "\033[100m";  // Grey background - empty/unplaced position
const string EMPTY_FG = "\033[90m";   // Grey foreground - empty position (grey bg + grey fg)

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
