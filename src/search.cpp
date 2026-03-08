#include <ankerl/unordered_dense.h>
#include <atomic>
#include <chrono>
#include <climits>
#include <deque>
#include <functional>
#include <iostream>
#include <mutex>
#include <robin_hood.h>
#include <thread>
#include <vector>
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include "config.h"

#ifdef _WIN32
// If system memory use is above 80%, trim our working set so the OS pages us
// out and avoids pushing the system into a crash.
static void trimWorkingSetIfHighMemory() {
  MEMORYSTATUSEX mem = {};
  mem.dwLength = sizeof(mem);
  if (!GlobalMemoryStatusEx(&mem))
    return;
  if (mem.dwMemoryLoad > 80)
    SetProcessWorkingSetSize(GetCurrentProcess(), (SIZE_T)-1, (SIZE_T)-1);
}
#endif
#include "display.h"
#include "grid_config.h"
#include "search.h"
#include "state_utils.h"
#include "transition.h"

using namespace std;

// Global file stream for baseline output (only used for reverse BFS section)
ofstream *g_baselineFile = nullptr;
streambuf *g_originalCoutBuf = nullptr;

// ── Packed-transition embedding ──────────────────────────────────────────────
// Level 5 uses at most 102 bits of the 128-bit state, so bits 120-127 (the top
// byte of uint128_simd::u.hi) are always zero in a canonical state.  We store
// the 8-bit Transition::packed there, turning every map<StateType,Transition>
// into a set<StateType> — eliminating the value slot and its alignment padding.
//
//   sizeof(pair<uint128_simd, uint8_t>) = 32 bytes (16-byte alignment + 15 pad)
//   sizeof(uint128_simd)               = 16 bytes   → 2× reduction per entry
//
// Encoding: top byte of hi = Transition::packed  (255 = sentinel / start node)
// The queues and all pure-state operations (stateGet, canonical, …) are safe:
//   • stateGet/stateSet only touch bits for positions 0-33 (max bit = 101)
//   • canonical()'s applyPerm initialises result=0 and only writes positions
//     0..numPos-1, so it naturally strips the packed byte → output is always a
//     pure canonical state ready for a fresh encode.

static constexpr uint64_t PACKED_HI_SHIFT = 56; // bit-offset within hi word
// Mask that keeps all bits except the top byte of hi (i.e. keeps bits 0-119).
static const StateType STATE_MASK =
    StateType(0xFFFFFFFFFFFFFFFFull, 0x00FFFFFFFFFFFFFFull);

inline StateType encodeTransition(StateType pureState, Transition tr) {
  return pureState | StateType(0, (uint64_t)tr.packed << PACKED_HI_SHIFT);
}
inline Transition decodeTransition(StateType encoded) {
  Transition t;
  t.packed = static_cast<uint8_t>(encoded.u.hi >> PACKED_HI_SHIFT);
  return t;
}

// ── Hash / equality that ignore the packed top byte ──────────────────────────
// Used for all visited/parent sets so that find(pureState) matches any entry
// with the same state bits regardless of what transition is encoded.

struct StateHash {
  size_t operator()(const StateType &x) const noexcept {
    return size_t(x.u.lo ^ (x.u.hi & 0x00FFFFFFFFFFFFFFull));
  }
};
struct StateEqual {
  bool operator()(const StateType &a, const StateType &b) const noexcept {
    return a.u.lo == b.u.lo &&
           (a.u.hi & 0x00FFFFFFFFFFFFFFull) == (b.u.hi & 0x00FFFFFFFFFFFFFFull);
  }
};

// Hash specialization for uint128_simd used by robin_hood containers.
// Ignores the top byte (packed transition) so that pure-state lookups work.
namespace robin_hood {
template <> struct hash<uint128_simd> {
  size_t operator()(const uint128_simd &x) const noexcept {
    return size_t(x.u.lo ^ (x.u.hi & 0x00FFFFFFFFFFFFFFull));
  }
};
} // namespace robin_hood

// ── Dynamic storage type
// ────────────────────────────────────────────────────── Levels 1+2 (numPos ≤
// 16): state ≤ 48 bits; packed fits in top byte of uint64_t (bits 56-63).
// Storing 8 bytes instead of 16 halves memory again. Levels 3-5 (numPos > 16):
// keep current StateType storage (bits 120-127).
//
// All BFS state operations (stateGet, canonical, …) always use StateType.
// Only the visited/parent SET keys are stored as the narrower uint64_t.
// Conversion boundary is at insert/find time via the helpers below.

// Defined after g_numPos (declared below in the globals section).
inline bool usesStorage64();

// Pure StateType → storage key (strips any packed bits, narrows if needed).
template <typename S> inline S toStorageKey(StateType s);
template <> inline uint64_t toStorageKey<uint64_t>(StateType s) {
  return s.u.lo;
}
template <> inline StateType toStorageKey<StateType>(StateType s) { return s; }

// Encode pure state + transition into a storage element.
template <typename S> inline S encodeForStorage(StateType state, Transition tr);
template <>
inline uint64_t encodeForStorage<uint64_t>(StateType state, Transition tr) {
  return state.u.lo | ((uint64_t)tr.packed << 56);
}
template <>
inline StateType encodeForStorage<StateType>(StateType state, Transition tr) {
  return encodeTransition(state, tr);
}

// Decode storage element → pure canonical StateType (strips packed bits).
template <typename S> inline StateType storageToState(S s);
template <> inline StateType storageToState<uint64_t>(uint64_t s) {
  return StateType(s & 0x00FFFFFFFFFFFFFFull, 0);
}
template <> inline StateType storageToState<StateType>(StateType s) {
  return s & STATE_MASK;
}

// Decode storage element → Transition.
template <typename S> inline Transition storageTrans(S s);
template <> inline Transition storageTrans<uint64_t>(uint64_t s) {
  Transition t;
  t.packed = static_cast<uint8_t>(s >> 56);
  return t;
}
template <> inline Transition storageTrans<StateType>(StateType s) {
  return decodeTransition(s);
}

// Per-storage hash and equality (ignore packed bits in either representation).
template <typename S> struct StorageHash;
template <> struct StorageHash<uint64_t> {
  size_t operator()(uint64_t x) const noexcept {
    return std::hash<uint64_t>{}(x & 0x00FFFFFFFFFFFFFFull);
  }
};
template <> struct StorageHash<StateType> : StateHash {};

template <typename S> struct StorageEqual;
template <> struct StorageEqual<uint64_t> {
  bool operator()(uint64_t a, uint64_t b) const noexcept {
    return (a & 0x00FFFFFFFFFFFFFFull) == (b & 0x00FFFFFFFFFFFFFFull);
  }
};
template <> struct StorageEqual<StateType> : StateEqual {};

// ── Hot-path lookup tables ──────────────────────────────────────────────────
// Flat fixed-size arrays eliminate vector pointer indirection in the inner
// loop. Max grid 7×7 = 49 positions, max 4 orthogonal neighbors each. Neighbor
// bit-shifts are precomputed (neighborIndex * 3) to avoid multiplies in
// canPlace. All positions are padded to exactly 4 entries; unused entries use
// shift 63 which reads bits 63-65 of the uint64_t — always 0 (type 0/empty),
// so they never set any bit that placement checks require (bits 1+).

static int g_numPos;
static int g_neighborShifts[49][4]; // precomputed: neighborIndex * 3
static int g_neighborCount[49];
static int g_maxType[49];
static StateType g_baseState; // all positions set to type 1 (Blue)
// Type T requires neighbors of types 1..T-1: mask = bits 1..T-1 set
static const uint32_t g_requiredMask[6] = {0, 0, 2, 6, 14, 30};

// usesStorage64() forward-declared above; defined here after g_numPos.
inline bool usesStorage64() { return g_numPos <= 16; }

// Standalone helpers (no lambdas) so profiler shows named functions.
static int countTowersLeft(StateType s) {
  int towersLeft = 0;
  int targetType = PATH_FROM_EMPTY ? 0 : 1;
  for (int pos = 0; pos < g_numPos; pos++) {
    if (stateGet(s, pos) != targetType)
      towersLeft++;
  }
  return towersLeft;
}

// ── Symmetry ────────────────────────────────────────────────────────────────
// Auto-detects valid symmetries by trying all 8 dihedral transforms against
// actual grid positions. Keeps only those where every position maps to another
// valid position. Identity is always first (index 0).

static int g_numSymm;
static int g_symmPerm[8][49]; // g_symmPerm[s][i] = image of position i under
                              // symmetry s

inline StateType applyPerm(StateType s, const int *perm) {
  StateType result = 0;
  for (int i = 0; i < g_numPos; i++)
    result |= StateType((s >> (i * 3)) & 7) << (perm[i] * 3);
  return result;
}

// Canonical form = minimum StateType value over all symmetry images.
// The minimum is a unique representative of the equivalence class.
inline StateType canonical(StateType s) {
  if (g_numSymm == 1)
    return s; // No symmetries, skip computation
  StateType best = s;
  for (int k = 1; k < g_numSymm; k++) {
    StateType t = applyPerm(s, g_symmPerm[k]);
    if (t < best) {
      best = t;
      if (best == g_baseState)
        return g_baseState; // Early exit: can't get smaller than all-blue
    }
  }
  return best;
}

// Returns {canonical(s), child_pos} where child_pos is where position `pos`
// ends up after the canonical permutation.  Used when building a Transition so
// we can store only (child_pos, prev_type) instead of the full parent state.
inline pair<StateType, int> canonicalAndChildPos(StateType s, int pos) {
  if (g_numSymm == 1)
    return {s, pos};
  StateType best = s;
  int bestK = 0;
  for (int k = 1; k < g_numSymm; k++) {
    StateType t = applyPerm(s, g_symmPerm[k]);
    if (t < best) {
      best = t;
      bestK = k;
      if (best == g_baseState)
        break;
    }
  }
  return {best, g_symmPerm[bestK][pos]};
}

static void buildSymmetries() {
  int R = GRID.rows, C = GRID.cols;
  // All 8 elements of the dihedral group D4:
  //   identity, 180° rot, 90° CW, 90° CCW, H-flip, V-flip, diagonal,
  //   anti-diagonal
  using T = function<pair<int, int>(int, int)>;
  T transforms[8] = {
      [](int r, int c) { return make_pair(r, c); },
      [R, C](int r, int c) { return make_pair(R - 1 - r, C - 1 - c); },
      [R](int r, int c) { return make_pair(c, R - 1 - r); },
      [C](int r, int c) { return make_pair(C - 1 - c, r); },
      [C](int r, int c) { return make_pair(r, C - 1 - c); },
      [R](int r, int c) { return make_pair(R - 1 - r, c); },
      [](int r, int c) { return make_pair(c, r); },
      [R, C](int r, int c) { return make_pair(C - 1 - c, R - 1 - r); },
  };

  g_numSymm = 0;
  for (int t = 0; t < 8; t++) {
    int perm[49];
    bool valid = true;
    for (int i = 0; i < g_numPos && valid; i++) {
      auto [r, c] = GRID.positions[i];
      auto [nr, nc] = transforms[t](r, c);
      auto it = GRID.coordToIdx.find({nr, nc});
      if (it == GRID.coordToIdx.end())
        valid = false;
      else
        perm[i] = it->second;
    }
    if (valid) {
      for (int i = 0; i < g_numPos; i++)
        g_symmPerm[g_numSymm][i] = perm[i];
      g_numSymm++;
    }
  }
  cout << "Symmetries detected: " << g_numSymm << "  (search space reduced ~"
       << g_numSymm << "x)\n";
}

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
  g_baseState = StateType(0);
  for (int i = 0; i < g_numPos; i++)
    g_baseState = stateSet(g_baseState, i, 1);
}

// ── Placement rule ───────────────────────────────────────────────────────────
// Tower of type T at pos requires one neighbor of each type 1..T-1.
// neighborPresent() computes the bitmask once per (state, pos); canPlace()
// then checks it cheaply for each candidate type.

inline uint32_t neighborPresent(StateType state, int pos) {
  uint32_t present = 0;
  const int *shifts = g_neighborShifts[pos];
  int count = g_neighborCount[pos];
  for (int i = 0; i < count; i++)
    present |= 1u << static_cast<int>((state >> shifts[i]) & 7);
  return present;
}

inline bool canPlace(uint32_t present, int type) {
  if (type <= 1)
    return true;
  return (present & g_requiredMask[type]) == g_requiredMask[type];
}

// ── Path reconstruction helpers ──────────────────────────────────────────────
static int scoreGrid(StateType state) {
  int score = 0;
  for (int i = 0; i < g_numPos; i++)
    score += stateGet(state, i) - 1;
  return score;
}

// expandPathForEmptyBoard: converts a path from the solver's "all-1 = blue"
// model to a path that starts from a truly empty board.
//
// For each step (pos, newType), before we can upgrade pos, we need:
//   1. pos itself to be explicitly placed (as blue) if not already
//   2. For each required neighbor type 1..newType-1, at least one neighbor
//      must be explicitly placed with that type.
//      - Neighbors with type > 1 were placed by earlier explicit steps.
//      - Neighbors with type == 1 are "implicit blue" and need an explicit
//        "place blue" step inserted here.
//
// Tracks a physical (non-canonical) state so neighbor lookups against
// GRID.neighbors[] remain consistent with position indices.
static vector<pair<int, int>>
expandPathForEmptyBoard(const vector<pair<int, int>> &steps) {
  vector<pair<int, int>> out;
  out.reserve(steps.size() * 2);
  uint64_t placed = 0; // bit i = physical position i has been explicitly placed

  // Track both states exactly like printExpandedPath does.
  StateType canonState = g_baseState;
  StateType physState = g_baseState;

  for (auto [pos, newType] : steps) {
    // Find which symmetry maps physState to its canonical form.
    int usedSymm = 0;
    for (int k = 0; k < g_numSymm; k++) {
      if (applyPerm(physState, g_symmPerm[k]) == canonState) {
        usedSymm = k;
        break;
      }
    }
    // perm[i]    = canonical position of physical position i
    // invPerm[j] = physical position of canonical position j
    int invPerm[49];
    for (int i = 0; i < g_numPos; i++)
      invPerm[g_symmPerm[usedSymm][i]] = i;
    int physPos = invPerm[pos];

    // For each required neighbor type req in 1..newType-1, ensure at least one
    // placed physical neighbor has that type. Unplaced neighbors with type 1
    // are "implicit blue"; insert an explicit place-blue step if needed.
    for (int req = 1; req < newType; req++) {
      bool satisfied = false;
      int candidatePhysNb = -1;
      for (int physNb : GRID.neighbors[physPos]) {
        int nbType = stateGet(physState, physNb);
        if (nbType == req) {
          if ((placed >> physNb) & 1) {
            satisfied = true;
            break;
          } else if (req == 1) {
            candidatePhysNb = physNb;
          }
        }
      }
      if (!satisfied && candidatePhysNb != -1) {
        // Express the placement in canonical space for printExpandedPath.
        // Blue placements don't change canonState (already 1), so usedSymm
        // stays stable for this step.
        int canonNb = g_symmPerm[usedSymm][candidatePhysNb];
        out.push_back({canonNb, 1});
        placed |= (1ull << candidatePhysNb);
      }
    }

    out.push_back({pos, newType});
    placed |= (1ull << physPos);
    canonState = canonical(stateSet(canonState, pos, newType));
    physState = stateSet(physState, physPos, newType);
  }

  return out;
}

static void advanceStateForDisplay(StateType &canonState,
                                   StateType &displayState,
                                   uint64_t &placedMask, int pos, int newType) {
  int usedSymm = 0;
  for (int k = 0; k < g_numSymm; k++) {
    if (applyPerm(displayState, g_symmPerm[k]) == canonState) {
      usedSymm = k;
      break;
    }
  }
  int inversePerm[49];
  for (int i = 0; i < g_numPos; i++)
    inversePerm[g_symmPerm[usedSymm][i]] = i;
  int physPos = inversePerm[pos];
  canonState = canonical(stateSet(canonState, pos, newType));
  displayState = stateSet(displayState, physPos, newType);
  placedMask |= (1ull << physPos);
}

static void printExpandedPath(const vector<pair<int, int>> &expandedSteps,
                              StateType startState) {
  const bool fromEmpty = (startState == StateType(0));
  StateType canonState = startState;
  StateType displayState = startState;
  uint64_t placedMask = 0;

  // Pre-pass to compute total steps (and trailing implicit blues when from
  // all-1).
  for (auto [pos, newType] : expandedSteps)
    advanceStateForDisplay(canonState, displayState, placedMask, pos, newType);

  int totalSteps = (int)expandedSteps.size();
  if (!fromEmpty) {
    int missingCount = 0;
    for (int i = 0; i < g_numPos; i++) {
      if (!((placedMask >> i) & 1))
        missingCount++;
    }
    totalSteps += missingCount;
  }

  // Reset for actual print pass.
  canonState = startState;
  displayState = startState;
  placedMask = 0;

  cout << "\n--- Build path (" << totalSteps << " steps) ---\n\n";
  cout << "Start (" << (fromEmpty ? "empty" : "all blue") << " board):\n";
  printGrid(startState, uint64_t(0));
  int validSteps = 0, invalidSteps = 0;
  int stepNum = 1;

  for (auto [pos, newType] : expandedSteps) {
    // pos is a position index in canonical space.
    // Find which symmetry maps our displayState to its canonical form,
    // then invert it to find the physical position.
    int usedSymm = 0;
    for (int k = 0; k < g_numSymm; k++) {
      if (applyPerm(displayState, g_symmPerm[k]) == canonState) {
        usedSymm = k;
        break;
      }
    }
    // Build inverse permutation: inversePerm[canonical_pos] = physical_pos
    int inversePerm[49];
    for (int i = 0; i < g_numPos; i++)
      inversePerm[g_symmPerm[usedSymm][i]] = i;

    int physPos = inversePerm[pos];
    bool wasPlaced = (placedMask >> physPos) & 1;
    int oldType = wasPlaced ? stateGet(displayState, physPos) : -1;

    // Validity check in canonical space
    uint32_t present = neighborPresent(canonState, pos);
    bool valid = canPlace(present, newType);
    if (valid)
      validSteps++;
    else
      invalidSteps++;

    advanceStateForDisplay(canonState, displayState, placedMask, pos, newType);

    auto [r, c] = GRID.positions[physPos];
    cout << "\nStep " << stepNum << ": pos " << physPos << " (" << r << "," << c
         << ")  ";
    if (oldType == -1) {
      cout << "empty -> " << COLORS[newType] << BLACK_TEXT
           << COLOR_NAMES[newType] << RESET << " (" << newType << ")";
    } else {
      cout << COLORS[oldType] << BLACK_TEXT << COLOR_NAMES[oldType] << RESET
           << " (" << oldType << ") -> " << COLORS[newType] << BLACK_TEXT
           << COLOR_NAMES[newType] << RESET << " (" << newType << ")";
    }
    cout << (valid ? "  \033[32m\u2713\033[0m"
                   : "  \033[31m\u2717 INVALID\033[0m")
         << "\n";
    printGrid(displayState, placedMask);
    stepNum++;
  }

  // When path started from all-blue: any position never explicitly placed is an
  // implicit blue; show a final step for each.
  if (!fromEmpty) {
    for (int i = 0; i < g_numPos; i++) {
      if (!((placedMask >> i) & 1)) {
        placedMask |= (1ull << i);
        auto [r, c] = GRID.positions[i];
        cout << "\nStep " << stepNum << ": pos " << i << " (" << r << "," << c
             << ")  ";
        cout << "empty -> " << COLORS[1] << BLACK_TEXT << COLOR_NAMES[1]
             << RESET << " (1)";
        cout << "  \033[32m\u2713\033[0m"
             << "\n";
        printGrid(displayState, placedMask);
        stepNum++;
        validSteps++;
      }
    }
  }

  int reachedScore = scoreGrid(displayState);
  cout << "\nBuild path: " << totalSteps << " steps, " << validSteps
       << " valid, " << invalidSteps << " invalid";
  if (invalidSteps == 0) {
    cout << " \033[32m\u2014 ALL VALID\033[0m";
  } else {
    cout << " \033[31m\u2014 HAS INVALID STEPS\033[0m";
  }
  cout << "\nReached score: " << reachedScore << "\n";
}

template <typename S>
inline void
printBuildPath(StateType bestState,
               const robin_hood::unordered_set<S, StorageHash<S>,
                                               StorageEqual<S>> &visited) {
  vector<pair<int, int>>
      steps; // (pos, newType) — pos in state *before* the move
  StateType cur = bestState;
  while (true) {
    auto tr = storageTrans<S>(*visited.find(toStorageKey<S>(cur)));
    if (tr.getPos() == -1)
      break;
    int childPos = tr.getPos();
    int newType = stateGet(cur, childPos);
    StateType parent = canonical(stateSet(cur, childPos, tr.getType()));
    // Replay expects (pos, newType) where pos is in the state before the move.
    // childPos is in cur (child); find pos in parent s.t. stateSet(parent, pos,
    // newType) -> cur.
    int posInParent = -1;
    for (int p = 0; p < g_numPos; p++) {
      if (stateGet(parent, p) == newType)
        continue;
      if (canonical(stateSet(parent, p, newType)) == cur) {
        posInParent = p;
        break;
      }
    }
    if (posInParent < 0)
      posInParent = childPos; // fallback if no symmetry (g_numSymm == 1)
    steps.push_back({posInParent, newType});
    cur = parent;
  }
  reverse(steps.begin(), steps.end());
  if (PATH_FROM_EMPTY)
    printExpandedPath(steps, StateType(0));
  else
    printExpandedPath(expandPathForEmptyBoard(steps), g_baseState);
}

struct SearchResult {
  StateType bestState;
  int bestScore;
  long long explored;
  size_t uniqueStates;
};

// BFS from base state (all blue) stopping as soon as target is reached.
// Returns a set of encoded states (transition packed into top byte) for path
// reconstruction.
inline robin_hood::unordered_set<StateType> bfsShortestPath(StateType target) {
  robin_hood::unordered_set<StateType> path;
  path.reserve(1ull << 24); // ~16M elements
  struct Entry {
    StateType state;
    int score;
  };
  vector<Entry> q;
  size_t head = 0;

  StateType bfsStart = PATH_FROM_EMPTY ? StateType(0) : g_baseState;
  path.emplace(encodeTransition(bfsStart, {})); // sentinel: packed=255
  q.push_back({bfsStart, 0});

  while (head < q.size()) {
    auto [cur, curScore] = q[head++];
    if (cur == target)
      break;
    for (int pos = 0; pos < g_numPos; pos++) {
      int curType = stateGet(cur, pos);
      uint32_t present = neighborPresent(cur, pos);
      for (int newType = 0; newType <= g_maxType[pos]; newType++) {
        if (!FORWARD_TRANSITIONS[curType][newType])
          continue;
        if (!canPlace(present, newType))
          continue;
        auto [next, childPos] =
            canonicalAndChildPos(stateSet(cur, pos, newType), pos);
        if (path.emplace(encodeTransition(next, Transition{(int8_t)childPos,
                                                           (int8_t)curType}))
                .second)
          q.push_back({next, curScore + (newType - curType)});
      }
    }
  }
  return path;
}

// Reverse BFS shared types and helpers (named functions for profiling).
using QueueEntry = std::pair<StateType, int>;

template <typename S> struct ReverseExpansion {
  S next;
  int nextTowers;
};

template <typename ParentSet>
static void pruneSideImpl(ParentSet &parentMap, std::deque<QueueEntry> &queue,
                          bool keepLowest) {
  using S = typename ParentSet::value_type;
  static constexpr size_t kPruneThreshold = 10000;
  static constexpr int kPruneLowDelta = 3;
  static constexpr int kPruneHighDelta = 3;
  if (parentMap.size() <= kPruneThreshold || queue.empty())
    return;

  int minTowers = g_numPos;
  int maxTowers = 0;
  for (const auto &entry : queue) {
    int tc = entry.second;
    minTowers = min(minTowers, tc);
    maxTowers = max(maxTowers, tc);
  }

  int cutoff =
      keepLowest ? (minTowers + kPruneLowDelta) : (maxTowers - kPruneHighDelta);
  vector<QueueEntry> toProcess;
  toProcess.reserve(queue.size());
  for (const auto &entry : queue) {
    if ((keepLowest && entry.second <= cutoff) ||
        (!keepLowest && entry.second >= cutoff)) {
      toProcess.push_back(entry);
    }
  }

  if (toProcess.empty()) {
    if (keepLowest)
      toProcess.push_back(
          *min_element(queue.begin(), queue.end(),
                       [](const QueueEntry &a, const QueueEntry &b) {
                         return a.second < b.second;
                       }));
    else
      toProcess.push_back(
          *max_element(queue.begin(), queue.end(),
                       [](const QueueEntry &a, const QueueEntry &b) {
                         return a.second > b.second;
                       }));
  }

  ParentSet pruned;
  pruned.reserve(toProcess.size() * 4);
  for (const auto &entry : toProcess)
    pruned.emplace(*parentMap.find(toStorageKey<S>(entry.first)));

  std::deque<QueueEntry> worklist(toProcess.begin(), toProcess.end());
  while (!worklist.empty()) {
    StateType cur = worklist.front().first;
    worklist.pop_front();
    auto it = parentMap.find(toStorageKey<S>(cur));
    if (it == parentMap.end())
      continue;
    auto tr = storageTrans<S>(*it);
    if (tr.getPos() == -1)
      continue;
    StateType parent = canonical(stateSet(cur, tr.getPos(), tr.getType()));
    if (parentMap.count(toStorageKey<S>(parent)) &&
        pruned.find(toStorageKey<S>(parent)) == pruned.end()) {
      pruned.emplace(*parentMap.find(toStorageKey<S>(parent)));
      worklist.push_back({parent, countTowersLeft(parent)});
    }
  }

  parentMap.swap(pruned);
  queue.assign(toProcess.begin(), toProcess.end());
}

template <typename S>
static void
reverseExpandChunkImpl(const QueueEntry *start, const QueueEntry *end,
                       vector<ReverseExpansion<S>> *out,
                       StateType pathStartState, int reverseBaseType) {
  const auto &backTable =
      PATH_FROM_EMPTY ? BACKWARD_TRANSITIONS_EMPTY : BACKWARD_TRANSITIONS;
  for (; start != end; ++start) {
    StateType cur = start->first;
    int curTowers = start->second;
    if (cur == pathStartState)
      continue;
    for (int pos = 0; pos < g_numPos; pos++) {
      int curType = stateGet(cur, pos);
      uint32_t present = neighborPresent(cur, pos);
      if (!canPlace(present, curType))
        continue;
      for (int prevType = 0; prevType <= g_maxType[pos]; prevType++) {
        if (!backTable[curType][prevType])
          continue;
        auto [next, childPos] =
            canonicalAndChildPos(stateSet(cur, pos, prevType), pos);
        int nextTowers = curTowers - (curType != reverseBaseType ? 1 : 0) +
                         (prevType != reverseBaseType ? 1 : 0);
        out->push_back({encodeForStorage<S>(next, Transition{(int8_t)childPos,
                                                             (int8_t)curType}),
                        nextTowers});
      }
    }
  }
}

template <typename S> struct ReverseTaskContext {
  std::deque<QueueEntry> *qReverse;
  ankerl::unordered_dense::set<S, StorageHash<S>, StorageEqual<S>>
      *parentReverse;
  StateType pathStartState;
  std::mutex *coutMutex;
  std::atomic<bool> *reverseFinished;
  size_t *totalReverseVisited;
  int reverseBaseType;
};

template <typename S> static void runReverseTask(ReverseTaskContext<S> *ctx) {
  std::deque<QueueEntry> &qReverse = *ctx->qReverse;
  auto &parentReverse = *ctx->parentReverse;
  StateType pathStartState = ctx->pathStartState;
  std::mutex &coutMutex = *ctx->coutMutex;
  std::atomic<bool> &reverseFinished = *ctx->reverseFinished;
  size_t &totalReverseVisited = *ctx->totalReverseVisited;
  int reverseBaseType = ctx->reverseBaseType;

  size_t depth = 0;
  while (!qReverse.empty()) {
    size_t layerCount = qReverse.size();
    {
      lock_guard<mutex> lock(coutMutex);
      cout << "[Reverse] Layer " << depth << ": expanding " << layerCount
           << " states ("
           << max(1, min(REVERSE_BFS_NUM_THREADS, (int)layerCount))
           << " threads)..." << endl;
    }
    auto layerStart = chrono::steady_clock::now();

    vector<QueueEntry> layer;
    layer.reserve(layerCount);
    for (size_t i = 0; i < layerCount; i++) {
      layer.push_back(qReverse.front());
      qReverse.pop_front();
    }

    int minTowersLeft = g_numPos;
    int maxTowersLeft = 0;
    bool reachedBase = false;
    for (const auto &entry : layer) {
      if (entry.first == pathStartState) {
        reachedBase = true;
        break;
      }
      minTowersLeft = min(minTowersLeft, entry.second);
      maxTowersLeft = max(maxTowersLeft, entry.second);
    }
    if (reachedBase) {
      qReverse.clear();
      auto layerEnd = chrono::steady_clock::now();
      double expandMs =
          chrono::duration_cast<chrono::milliseconds>(layerEnd - layerStart)
              .count();
      lock_guard<mutex> lock(coutMutex);
      cout << "[Reverse] Layer " << depth << " done in " << expandMs
           << " ms (expand " << expandMs << " ms). Reached base state." << endl;
      break;
    }

    const int numThreads =
        max(1, min(REVERSE_BFS_NUM_THREADS, (int)layer.size()));
    const size_t chunkSize =
        (layer.size() + (size_t)numThreads - 1) / (size_t)numThreads;
    vector<vector<ReverseExpansion<S>>> results((size_t)numThreads);
    for (int c = 0; c < numThreads; c++)
      results[c].reserve(layer.size() * 2);

    if (numThreads == 1) {
      reverseExpandChunkImpl<S>(layer.data(), layer.data() + layer.size(),
                                &results[0], pathStartState, reverseBaseType);
    } else {
      vector<thread> revThreads;
      for (int c = 0; c < numThreads; c++) {
        const QueueEntry *chunkStart = layer.data() + (size_t)c * chunkSize;
        const QueueEntry *chunkEnd =
            layer.data() + min((size_t)(c + 1) * chunkSize, layer.size());
        if (chunkStart < chunkEnd)
          revThreads.emplace_back(reverseExpandChunkImpl<S>, chunkStart,
                                  chunkEnd, &results[c], pathStartState,
                                  reverseBaseType);
      }
      for (auto &t : revThreads)
        t.join();
    }

    {
      size_t inserted = 0;
      constexpr size_t trimCheckInterval = 50000;
      for (int c = 0; c < numThreads; c++) {
        for (const auto &r : results[c]) {
          if (parentReverse.emplace(r.next).second) {
            totalReverseVisited++;
            qReverse.push_back({storageToState<S>(r.next), r.nextTowers});
            if (++inserted % trimCheckInterval == 0) {
#ifdef _WIN32
              trimWorkingSetIfHighMemory();
#endif
            }
          }
        }
      }
    }

    auto layerEnd = chrono::steady_clock::now();
    double expandMs =
        chrono::duration_cast<chrono::milliseconds>(layerEnd - layerStart)
            .count();
    double pruneMs = 0.0;
    if (!qReverse.empty()) {
      auto pruneStart = chrono::steady_clock::now();
      pruneSideImpl(parentReverse, qReverse, true);
      pruneMs = chrono::duration_cast<chrono::milliseconds>(
                    chrono::steady_clock::now() - pruneStart)
                    .count();
    }
    double totalMs = expandMs + pruneMs;
    {
      lock_guard<mutex> lock(coutMutex);
      cout << "[Reverse] Layer " << depth << " done in " << totalMs
           << " ms (expand " << expandMs << " ms, prune " << pruneMs
           << " ms). Queue size now " << qReverse.size() << ", map size "
           << parentReverse.size() << "." << endl;
    }
#ifdef _WIN32
    trimWorkingSetIfHighMemory();
#endif
    depth++;
  }
  {
    lock_guard<mutex> lock(coutMutex);
    cout << "[Reverse] Queue empty. Reverse BFS finished." << endl;
  }
  reverseFinished = true;
}

template <typename S> struct ForwardTaskContext {
  std::deque<QueueEntry> *qForward;
  ankerl::unordered_dense::set<S, StorageHash<S>, StorageEqual<S>>
      *parentForward;
  StateType pathStartState;
  int reverseBaseType;
  std::atomic<bool> *reverseFinished;
  std::mutex *coutMutex;
};

template <typename S> static void runForwardTask(ForwardTaskContext<S> *ctx) {
  std::deque<QueueEntry> &qForward = *ctx->qForward;
  auto &parentForward = *ctx->parentForward;
  int reverseBaseType = ctx->reverseBaseType;
  std::atomic<bool> &reverseFinished = *ctx->reverseFinished;
  std::mutex &coutMutex = *ctx->coutMutex;

  size_t forwardDepth = 0;
  while (!qForward.empty()) {
    if (reverseFinished) {
      lock_guard<mutex> lock(coutMutex);
      cout << "[Forward] Reverse finished; stopping immediately." << endl;
      break;
    }
    size_t forwardLayerCount = qForward.size();
    {
      lock_guard<mutex> lock(coutMutex);
      cout << "[Forward] Layer " << forwardDepth << ": expanding "
           << forwardLayerCount << " states..." << endl;
    }
    auto forwardLayerStart = chrono::steady_clock::now();
    int minForwardTowersLeft = g_numPos;
    int maxForwardTowersLeft = 0;
    bool anyForwardProcessed = false;
    bool stoppedByReverse = false;

    for (size_t i = 0; i < forwardLayerCount; i++) {
      if (reverseFinished) {
        stoppedByReverse = true;
        break;
      }
      StateType cur = qForward.front().first;
      int curTowers = qForward.front().second;
      qForward.pop_front();
      if (!anyForwardProcessed) {
        minForwardTowersLeft = curTowers;
        maxForwardTowersLeft = curTowers;
        anyForwardProcessed = true;
      } else {
        minForwardTowersLeft = min(minForwardTowersLeft, curTowers);
        maxForwardTowersLeft = max(maxForwardTowersLeft, curTowers);
      }

      for (int pos = 0; pos < g_numPos; pos++) {
        int curType = stateGet(cur, pos);
        uint32_t present = neighborPresent(cur, pos);
        for (int newType = 0; newType <= g_maxType[pos]; newType++) {
          if (!FORWARD_TRANSITIONS[curType][newType])
            continue;
          if (!canPlace(present, newType))
            continue;
          auto [next, childPos] =
              canonicalAndChildPos(stateSet(cur, pos, newType), pos);
          if (parentForward
                  .emplace(encodeForStorage<S>(
                      next, Transition{(int8_t)childPos, (int8_t)curType}))
                  .second) {
            int nextTowers = curTowers - (curType != reverseBaseType ? 1 : 0) +
                             (newType != reverseBaseType ? 1 : 0);
            qForward.push_back({storageToState<S>(next), nextTowers});
          }
        }
      }
    }

    auto forwardLayerEnd = chrono::steady_clock::now();
    double expandMs = chrono::duration_cast<chrono::milliseconds>(
                          forwardLayerEnd - forwardLayerStart)
                          .count();
    double pruneMs = 0.0;
    if (!stoppedByReverse && !qForward.empty()) {
      auto pruneStart = chrono::steady_clock::now();
      pruneSideImpl(parentForward, qForward, false);
      pruneMs = chrono::duration_cast<chrono::milliseconds>(
                    chrono::steady_clock::now() - pruneStart)
                    .count();
    }
    double totalMs = expandMs + pruneMs;
    {
      lock_guard<mutex> lock(coutMutex);
      cout << "[Forward] Layer " << forwardDepth << " done in " << totalMs
           << " ms (expand " << expandMs << " ms, prune " << pruneMs
           << " ms). Queue size now " << qForward.size() << ", map size "
           << parentForward.size();
      if (stoppedByReverse)
        cout << " (stopped mid-layer: reverse finished)";
      cout << "." << endl;
    }
    if (stoppedByReverse)
      break;
    forwardDepth++;
  }
  {
    lock_guard<mutex> lock(coutMutex);
    cout << "[Forward] Forward BFS finished (queue empty or reverse finished)."
         << endl;
  }
}

template <typename S>
static int forwardDepthOfImpl(
    const ankerl::unordered_dense::set<S, StorageHash<S>, StorageEqual<S>>
        &parentForward,
    StateType s) {
  int d = 0;
  StateType cur = s;
  while (parentForward.count(toStorageKey<S>(cur))) {
    auto tr = storageTrans<S>(*parentForward.find(toStorageKey<S>(cur)));
    if (tr.getPos() == -1)
      break;
    cur = canonical(stateSet(cur, tr.getPos(), tr.getType()));
    d++;
  }
  return d;
}

template <typename S>
static int reverseDepthOfImpl(
    const ankerl::unordered_dense::set<S, StorageHash<S>, StorageEqual<S>>
        &parentReverse,
    StateType s) {
  int d = 0;
  StateType cur = s;
  while (parentReverse.count(toStorageKey<S>(cur))) {
    auto tr = storageTrans<S>(*parentReverse.find(toStorageKey<S>(cur)));
    if (tr.getPos() == -1)
      break;
    cur = canonical(stateSet(cur, tr.getPos(), tr.getType()));
    d++;
  }
  return d;
}

// Reverse BFS from a final state back to the path start (empty or all-blue
// depending on PATH_FROM_EMPTY). Used by mode 1 (Sliding Window) for path
// reconstruction since the search itself doesn't store transitions.
// S = uint64_t for levels 1+2 (state ≤ 48 bits), StateType for 3-5.
template <typename S>
static vector<pair<int, int>>
buildPathFromReverseBfsImpl(StateType finalState) {
  const StateType pathStartState = PATH_FROM_EMPTY ? StateType(0) : g_baseState;
  const int reverseBaseType = PATH_FROM_EMPTY ? 0 : 1;

#ifdef _WIN32
  // Hint that our working set can be trimmed so the OS pages us out earlier
  // instead of keeping RAM at 99% and only paging under crisis.
  SetProcessWorkingSetSize(GetCurrentProcess(), (SIZE_T)-1, (SIZE_T)-1);
#endif

  ankerl::unordered_dense::set<S, StorageHash<S>, StorageEqual<S>>
      parentReverse;
  parentReverse.reserve(1ull << 26); // ~64M elements
  ankerl::unordered_dense::set<S, StorageHash<S>, StorageEqual<S>>
      parentForward;
  if (RUN_PARALLEL_FORWARD_BFS)
    parentForward.reserve(1ull << 26); // ~64M elements

  std::deque<QueueEntry> qReverse;
  std::deque<QueueEntry> qForward;
  if (RUN_PARALLEL_FORWARD_BFS) {
    parentForward.emplace(encodeForStorage<S>(pathStartState, {}));
    qForward.push_back({pathStartState, countTowersLeft(pathStartState)});
  }
  parentReverse.emplace(encodeForStorage<S>(finalState, {}));
  qReverse.push_back({finalState, countTowersLeft(finalState)});

  std::mutex coutMutex;
  std::atomic<bool> reverseFinished{false};
  size_t totalReverseVisited =
      0; // total states ever discovered by reverse (before pruning)

  {
    lock_guard<mutex> lock(coutMutex);
    if (RUN_PARALLEL_FORWARD_BFS)
      cout << "[Main] Starting reverse BFS (to completion) and forward BFS in "
              "parallel. Forward will stop as soon as reverse finishes."
           << endl;
    else
      cout << "[Main] Reverse BFS only (parallel forward disabled)." << endl;
  }

  ReverseTaskContext<S> revCtx{
      &qReverse,        &parentReverse,       pathStartState, &coutMutex,
      &reverseFinished, &totalReverseVisited, reverseBaseType};
  thread tReverse(runReverseTask<S>, &revCtx);
  thread tForward;
  std::optional<ForwardTaskContext<S>> fwdCtx;
  if (RUN_PARALLEL_FORWARD_BFS) {
    fwdCtx =
        ForwardTaskContext<S>{&qForward,       &parentForward,   pathStartState,
                              reverseBaseType, &reverseFinished, &coutMutex};
    tForward = thread(runForwardTask<S>, &*fwdCtx);
  }
  tReverse.join();
  if (RUN_PARALLEL_FORWARD_BFS) {
    {
      lock_guard<mutex> lock(coutMutex);
      cout << "[Main] Reverse finished. Not waiting for forward to complete; "
              "stopping forward and finding best intercept with current state."
           << endl;
    }
    tForward.join();
  }

  // Redirect output to file if baseline mode is active (after BFS table
  // completes)
  if (g_baselineFile && g_baselineFile->is_open()) {
    g_originalCoutBuf = cout.rdbuf();
    cout.rdbuf(g_baselineFile->rdbuf());
  }

  cout << "\nReverse BFS visited " << totalReverseVisited
       << " states (map size after pruning: " << parentReverse.size() << ").\n";

  vector<pair<int, int>> steps;

  if (!RUN_PARALLEL_FORWARD_BFS) {
    // Walk from pathStartState toward finalState following the reverse BFS
    // tree.
    if (parentReverse.count(toStorageKey<S>(pathStartState)) == 0) {
      lock_guard<mutex> lock(coutMutex);
      cout << "[Main] Reverse BFS did not reach path start. Cannot build path."
           << endl;
      return {};
    }
    StateType cur = pathStartState;
    while (cur != finalState) {
      auto tr = storageTrans<S>(*parentReverse.find(toStorageKey<S>(cur)));
      StateType parent = canonical(stateSet(cur, tr.getPos(), tr.getType()));
      bool foundStep = false;
      int stepTypeMin = PATH_FROM_EMPTY ? 0 : 1;
      for (int p = 0; p < g_numPos && !foundStep; p++) {
        for (int t = stepTypeMin; t <= g_maxType[p] && !foundStep; t++) {
          if (t == stateGet(cur, p))
            continue;
          if (canonical(stateSet(cur, p, t)) == parent) {
            steps.push_back({p, t});
            foundStep = true;
          }
        }
      }
      cur = parent;
    }
    return steps;
  }

  // Meet-in-the-middle: find best intercept and build path from both sets.
  StateType meetState = pathStartState;
  int bestTotal = INT_MAX;
  size_t interceptCount = 0;
  for (const auto &elem : parentForward) {
    StateType s = storageToState<S>(elem);
    if (parentReverse.count(toStorageKey<S>(s)) == 0)
      continue;
    interceptCount++;
    int fd = forwardDepthOfImpl<S>(parentForward, s);
    int rd = reverseDepthOfImpl<S>(parentReverse, s);
    int total = fd + rd;
    if (total < bestTotal) {
      bestTotal = total;
      meetState = s;
    }
  }

  {
    lock_guard<mutex> lock(coutMutex);
    cout << "[Main] Found " << interceptCount
         << " states in both maps. Best intercept has forward depth "
         << forwardDepthOfImpl<S>(parentForward, meetState)
         << ", reverse depth "
         << reverseDepthOfImpl<S>(parentReverse, meetState) << ", total "
         << bestTotal << " steps." << endl;
  }

  if (interceptCount == 0) {
    lock_guard<mutex> lock(coutMutex);
    cout << "[Main] No intercept: no state reached by both searches. Cannot "
            "reconstruct path."
         << endl;
    return {};
  }

  vector<pair<int, int>> stepsForward;
  {
    StateType cur = meetState;
    while (true) {
      auto tr = storageTrans<S>(*parentForward.find(toStorageKey<S>(cur)));
      if (tr.getPos() == -1)
        break;
      int childPos = tr.getPos();
      int newType = stateGet(cur, childPos);
      StateType parent = canonical(stateSet(cur, childPos, tr.getType()));
      int posInParent = -1;
      for (int p = 0; p < g_numPos; p++) {
        if (stateGet(parent, p) == newType)
          continue;
        if (canonical(stateSet(parent, p, newType)) == cur) {
          posInParent = p;
          break;
        }
      }
      if (posInParent < 0)
        posInParent = childPos;
      stepsForward.push_back({posInParent, newType});
      cur = parent;
    }
  }
  reverse(stepsForward.begin(), stepsForward.end());

  {
    StateType cur = meetState;
    while (true) {
      auto tr = storageTrans<S>(*parentReverse.find(toStorageKey<S>(cur)));
      if (tr.getPos() == -1)
        break;
      StateType parent = canonical(stateSet(cur, tr.getPos(), tr.getType()));
      bool foundStep = false;
      int stepTypeMin = PATH_FROM_EMPTY ? 0 : 1;
      for (int p = 0; p < g_numPos && !foundStep; p++) {
        for (int t = stepTypeMin; t <= g_maxType[p] && !foundStep; t++) {
          if (t == stateGet(cur, p))
            continue;
          if (canonical(stateSet(cur, p, t)) == parent) {
            steps.push_back({p, t});
            foundStep = true;
          }
        }
      }
      cur = parent;
    }
  }

  steps.insert(steps.begin(), stepsForward.begin(), stepsForward.end());
  return steps;
}

// Runtime dispatch: use uint64_t storage for levels 1+2, StateType for 3-5.
static vector<pair<int, int>> buildPathFromReverseBfs(StateType finalState) {
  if (usesStorage64())
    return buildPathFromReverseBfsImpl<uint64_t>(finalState);
  else
    return buildPathFromReverseBfsImpl<StateType>(finalState);
}

// ── Sliding Window Optimization ─────────────────────────────────────────────-
// Processes grid using 3x3 windows with stride=1 along a JPEG-like zigzag path
// Each window respects constraints from the entire grid (in+out of window)

// Build zigzag path starting from a given corner
// corner: 0=top-left, 1=top-right, 2=bottom-left, 3=bottom-right
static vector<int> buildZigzagPath(int corner = 0) {
  int R = GRID.rows, C = GRID.cols;
  vector<vector<pair<int, int>>> diagonals;

  // Group by diagonal - direction depends on corner
  for (int r = 0; r < R; r++) {
    for (int c = 0; c < C; c++) {
      pair<int, int> coord = {r, c};
      if (!GRID.coordToIdx.count(coord))
        continue;

      int diag;
      switch (corner) {
      case 0:
        diag = r + c;
        break; // top-left
      case 1:
        diag = r + (C - 1 - c);
        break; // top-right
      case 2:
        diag = (R - 1 - r) + c;
        break; // bottom-left
      case 3:
        diag = (R - 1 - r) + (C - 1 - c);
        break; // bottom-right
      default:
        diag = r + c;
        break;
      }
      if (diag >= (int)diagonals.size())
        diagonals.resize(diag + 1);
      diagonals[diag].push_back({r, c});
    }
  }

  vector<int> path;
  for (int d = 0; d < (int)diagonals.size(); d++) {
    auto &diag = diagonals[d];
    // Alternate direction for each diagonal (zigzag effect)
    if (d % 2 == 1)
      reverse(diag.begin(), diag.end());
    for (auto [r, c] : diag) {
      path.push_back(GRID.coordToIdx[{r, c}]);
    }
  }
  return path;
}

// BFS within a window: exhaustively explore ALL reachable states
// by modifying only positions within the window
static void exploreWindow(StateType startState, int windowCenterIdx,
                          vector<StateType> &newStates, int &bestScore,
                          StateType &bestState, long long &explored) {
  // Get kernel pattern window around position
  auto [cr, cc] = GRID.positions[windowCenterIdx];

  // Parse kernel pattern to get offsets
  int patternRows = (int)WINDOW_KERNEL_PATTERN.size();
  int patternCols = 0;
  for (const auto &row : WINDOW_KERNEL_PATTERN) {
    patternCols = max(patternCols, (int)row.size());
  }
  int centerR = patternRows / 2;
  int centerC = patternCols / 2;

  vector<int> windowPositions;
  for (int pr = 0; pr < patternRows; pr++) {
    for (int pc = 0; pc < (int)WINDOW_KERNEL_PATTERN[pr].size(); pc++) {
      if (WINDOW_KERNEL_PATTERN[pr][pc] == '.') {
        int dr = pr - centerR;
        int dc = pc - centerC;
        auto it = GRID.coordToIdx.find({cr + dr, cc + dc});
        if (it != GRID.coordToIdx.end()) {
          windowPositions.push_back(it->second);
        }
      }
    }
  }

  // Full BFS within window - explore ALL states reachable by modifying window
  // positions
  robin_hood::unordered_set<StateType> localVisited;
  localVisited.reserve(1ull << 24); // ~16M elements per window
  struct Entry {
    StateType state;
    int score;
    int depth;
  };
  vector<Entry> q;

  localVisited.insert(startState);
  q.push_back({startState, scoreGrid(startState), 0});
  size_t head = 0;

  while (head < q.size()) {
    auto [cur, curScore, depth] = q[head++];
    explored++;

    // Stop expanding this branch if we've reached max depth
    if (depth >= MAX_DEPTH)
      continue;

    if (curScore > bestScore) {
      bestScore = curScore;
      bestState = cur;
      cout << "New best: " << bestScore << "\n";
      printGrid(bestState);
    }

    // Try modifying each position in the window
    for (int pos : windowPositions) {
      if (pos >= g_numPos)
        continue;
      int curType = stateGet(cur, pos);

      uint32_t present = neighborPresent(cur, pos);
      for (int newType = 0; newType <= g_maxType[pos]; newType++) {
        if (!FORWARD_TRANSITIONS_SLIDING_WINDOW[curType][newType])
          continue;
        if (!canPlace(present, newType))
          continue;

        StateType next = stateSet(cur, pos, newType);

        // Only track within this window to prevent infinite loops
        if (localVisited.insert(next).second) {
          int nextScore = curScore + (newType - curType);
          q.push_back({next, nextScore, depth + 1});
          newStates.push_back(next);
        }
      }
    }
  }
  // localVisited is destroyed here - no global state persists
}

static SearchResult runSlidingWindow() {
  SearchResult result{g_baseState, 0, 0, 0};
  result.bestState = g_baseState;
  result.bestScore = 0;
  result.explored = 0;

  cout << "Sliding window - sequential greedy approach:\n";
  for (const auto &row : WINDOW_KERNEL_PATTERN) {
    cout << "  " << row << "\n";
  }

  const char *cornerNames[] = {"top-left", "top-right", "bottom-left",
                               "bottom-right"};

  struct WindowStats {
    int posIdx;
    int row, col;
    int kernelSize;
    long long statesExplored;
    long long timeMs;
    bool improved;
    int scoreAfter;
  };

  StateType currentState = g_baseState;
  int currentScore = 0;
  int roundNum = 0;

  while (true) {
    roundNum++;
    int scoreAtRoundStart = currentScore;

    cout << "============ Round " << roundNum
         << " (starting score: " << currentScore << ") ============\n\n";

    int corner = 1; // top-right
    vector<int> zigzag = buildZigzagPath(corner);

    int bestScore = currentScore;
    StateType bestState = currentState;
    long long explored = 0;

    cout << "=== Pass " << (char)('A' + corner) << ": " << cornerNames[corner]
         << " (path length: " << zigzag.size() << ") ===\n";

    auto passStart = chrono::steady_clock::now();

    vector<WindowStats> stats;

    for (int windowIdx = 0; windowIdx < (int)zigzag.size(); windowIdx++) {
      vector<StateType> exploredStates;
      long long exploredBefore = explored;

      int posIdx = zigzag[windowIdx];
      auto [cr, cc] = GRID.positions[posIdx];

      auto t0 = chrono::steady_clock::now();

      exploreWindow(currentState, posIdx, exploredStates, bestScore, bestState,
                    explored);

      auto t1 = chrono::steady_clock::now();
      auto ms = chrono::duration_cast<chrono::milliseconds>(t1 - t0).count();

      bool improved = (bestState != currentState);
      if (improved) {
        currentState = bestState;
      }

      stats.push_back({posIdx, cr, cc, (int)exploredStates.size(),
                       explored - exploredBefore, ms, improved, bestScore});
    }

    currentScore = bestScore;

    auto passEnd = chrono::steady_clock::now();
    auto passMs =
        chrono::duration_cast<chrono::milliseconds>(passEnd - passStart)
            .count();

    // Print summary table
    cout << "\n┌─────┬─────────┬──────────────┬────────────────┬──────────┬────"
            "──┬───────"
            "┐\n";
    cout << "│ Pos │  (r,c)  │    Kernel    │     States     │   Time   │ Imp  "
            "│ "
            "Score │\n";
    cout << "├─────┼─────────┼──────────────┼────────────────┼──────────┼──────"
            "┼──────"
            "─┤\n";
    for (const auto &s : stats) {
      printf("│ %3d │ (%d, %d)  │ %12d │ %14lld │ %6lldms │  %s  │ %5d │\n",
             s.posIdx, s.row, s.col, s.kernelSize, s.statesExplored, s.timeMs,
             s.improved ? "Y" : "N", s.scoreAfter);
    }
    cout << "└─────┴─────────┴──────────────┴────────────────┴──────────┴──────"
            "┴──────"
            "─┘\n";
    cout << "Pass " << (char)('A' + corner) << " (" << cornerNames[corner]
         << "): score " << currentScore << " (" << passMs << "ms, " << explored
         << " states)\n\n";

    result.explored += explored;

    cout << "Round " << roundNum << " complete: score " << scoreAtRoundStart
         << " -> " << currentScore << "\n";

    if (currentScore <= scoreAtRoundStart) {
      cout << "No improvement in round " << roundNum << ", stopping.\n\n";
      break;
    } else {
      cout << "Improved by " << (currentScore - scoreAtRoundStart)
           << ", starting another round.\n\n";
    }
  }

  result.bestScore = currentScore;
  result.bestState = currentState;

  cout << "\nFinal score after " << roundNum
       << " round(s): " << result.bestScore << "\n";

  // ── Build path reconstruction: reverse BFS ─────────────────────────────
  StateType finalState = canonical(result.bestState);
  vector<pair<int, int>> buildPath = buildPathFromReverseBfs(finalState);
  if (buildPath.empty()) {
    cout << "\nBuild path not found by reverse BFS.\n";
    cout << "\nFinal score: " << result.bestScore << "\n";
    result.uniqueStates = 0;
    return result;
  }

  if (PATH_FROM_EMPTY)
    printExpandedPath(buildPath, StateType(0));
  else
    printExpandedPath(expandPathForEmptyBoard(buildPath), g_baseState);

  cout << "\nFinal score: " << result.bestScore << "\n";
  result.uniqueStates = 0;

  // Restore cout if it was redirected
  if (g_originalCoutBuf) {
    cout.rdbuf(g_originalCoutBuf);
    g_originalCoutBuf = nullptr;
  }

  return result;
}

// BFS over all reachable states.  S = uint64_t for levels 1+2, StateType for
// levels 3+.  Populates result and writes path reconstruction to cout (which
// may be redirected to the baseline file by the caller).
template <typename S> static void runBFSSearch(SearchResult &result) {
  robin_hood::unordered_set<S, StorageHash<S>, StorageEqual<S>> visited;
  visited.reserve(1ull << 24);

  struct Entry {
    StateType state;
    int score;
  };
  vector<Entry> q;
  size_t head = 0;

  StateType startState = PATH_FROM_EMPTY ? StateType(0) : g_baseState;
  visited.emplace(encodeForStorage<S>(startState, {})); // sentinel: no parent
  q.push_back({startState, 0});

  result.bestState = startState;
  result.bestScore = 0;
  result.explored = 0;

  while (head < q.size()) {
    auto [cur, curScore] = q[head++];
    result.explored++;

    if (curScore > result.bestScore) {
      result.bestScore = curScore;
      result.bestState = cur;
      cout << "New best score: " << result.bestScore << "\n";
      printGrid(cur);
    }

    for (int pos = 0; pos < g_numPos; pos++) {
      int curType = stateGet(cur, pos);
      uint32_t present = neighborPresent(cur, pos);
      for (int newType = 0; newType <= g_maxType[pos]; newType++) {
        if (!FORWARD_TRANSITIONS[curType][newType])
          continue;
        if (!canPlace(present, newType))
          continue;

        auto [next, childPos] =
            canonicalAndChildPos(stateSet(cur, pos, newType), pos);
        int nextScore = curScore + (newType - curType);
        if (visited
                .emplace(encodeForStorage<S>(
                    next, Transition{(int8_t)childPos, (int8_t)curType}))
                .second)
          q.push_back({next, nextScore});
      }
    }

    if (result.explored % 2000000 == 0)
      cout << "Explored: " << result.explored
           << " | Visited: " << visited.size()
           << " | Best: " << result.bestScore << "\r" << flush;
  }

  result.uniqueStates = visited.size();

  if (g_baselineFile && g_baselineFile->is_open()) {
    g_originalCoutBuf = cout.rdbuf();
    cout.rdbuf(g_baselineFile->rdbuf());
  }
  printBuildPath<S>(result.bestState, visited);
  if (g_originalCoutBuf) {
    cout.rdbuf(g_originalCoutBuf);
    g_originalCoutBuf = nullptr;
  }
}

// ── Public entrypoint for a single level/unlocked-type configuration ────────

void runSingleLevel() {
  cout << "Tower Blox Deluxe Solver\n";
  cout << "========================\n";
  cout << "Level: " << getCurrentLevel() << "\n";
  cout << "Tower types unlocked (1-" << MAX_UNLOCKED_TYPE << "):\n";
  for (int i = 1; i <= MAX_UNLOCKED_TYPE; i++)
    cout << "  " << COLORS[i] << BLACK_TEXT << i << RESET << " ("
         << COLOR_NAMES[i] << ")\n";
  if (SEARCH_MODE == 0)
    cout << "\nRunning BFS over all reachable states...\n\n";
  else
    cout << "\nRunning sliding window optimization (sequential greedy)...\n\n";

  buildLookupTables();
  buildSymmetries();

  SearchResult result{g_baseState, 0, 0, 0};

  if (SEARCH_MODE == 1) {
    result = runSlidingWindow();
  } else {
    // Dispatch to the appropriately-typed BFS (uint64_t for levels 1+2,
    // StateType for 3+).
    if (usesStorage64())
      runBFSSearch<uint64_t>(result);
    else
      runBFSSearch<StateType>(result);
  }

  cout << (SEARCH_MODE == 1 ? "\nSliding window complete.\n"
                            : "\nBFS complete.\n");
  cout << "States explored : " << result.explored << "\n";
  cout << "Unique states   : " << result.uniqueStates << "\n";
  cout << "Best score      : " << result.bestScore << "\n\n";
  cout << "Optimal layout:\n";
  printGrid(result.bestState);

  // Mode 0: path already printed by runBFSSearch above.
  // Mode 1: path already printed by runSlidingWindow.
}
