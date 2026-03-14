#pragma once

// Shared state, storage, and canonicalization used by search.cpp and
// search_bfs.cpp. Globals are declared here and defined in search.cpp.

#include "config.h"
#include "transition.h"

#include <ankerl/unordered_dense.h>

#include <cstddef>
#include <cstdint>
#include <utility>

// ── Globals (defined in search.cpp) ─────────────────────────────────────────
extern int g_numPos;
extern int g_neighborShifts[49][4];
extern int g_neighborCount[49];
extern int g_maxType[49];
extern StateType g_baseState;
extern const uint32_t g_requiredMask[6];
extern int g_numSymm;
extern int g_symmPerm[8][49];
extern int g_symmPermInv[8][49];

// ════════════════════════════════════════════════════════════════════════════
// State packing utilities
// ════════════════════════════════════════════════════════════════════════════
// Each position stores a value 0-5 (0=empty, 1-5=tower types, fits in 3 bits).
// Up to 21 positions × 3 bits = 63 bits → fits in uint64_t.

// Precomputed masks (defined in search.cpp, sized for 49 max positions):
// clearMask[pos]  = all 1s except the 3 bits for pos are 0.
// setMask[pos][v] = only the 3 bits for pos set to v (0..5).
extern StateType g_stateClearMask[49];
extern StateType g_stateSetMask[49][6];

inline StateType stateSet(StateType s, int pos, int val) {
  return (s & g_stateClearMask[pos]) | g_stateSetMask[pos][val];
}

inline int stateGet(StateType s, int pos) { return (s >> (pos * 3)) & 7; }

// ── Packed-transition embedding ─────────────────────────────────────────────
static constexpr uint64_t PACKED_HI_SHIFT = 56;
static const StateType STATE_MASK = StateType(0xFFFFFFFFFFFFFFFFull, 0x00FFFFFFFFFFFFFFull);

inline StateType encodeTransition(StateType pureState, Transition tr) {
  return pureState | StateType(0, (uint64_t)tr.packed << PACKED_HI_SHIFT);
}
inline Transition decodeTransition(StateType encoded) {
  Transition t;
  t.packed = static_cast<uint8_t>(encoded.u.hi >> PACKED_HI_SHIFT);
  return t;
}

// ── Hash / equality (ignore packed top byte) ───────────────────────────────
// High-quality hash (is_avalanching) so ankerl skips extra mixing; wyhash for
// good distribution.
struct StateHash {
  using is_avalanching = void;
  size_t operator()(const StateType &x) const noexcept {
    uint64_t parts[2] = {x.u.lo, x.u.hi & 0x00FFFFFFFFFFFFFFULL};
    return size_t(ankerl::unordered_dense::detail::wyhash::hash(parts, 16));
  }
};
struct StateEqual {
  bool operator()(const StateType &a, const StateType &b) const noexcept {
    return a.u.lo == b.u.lo && (a.u.hi & 0x00FFFFFFFFFFFFFFull) == (b.u.hi & 0x00FFFFFFFFFFFFFFull);
  }
};

// StateHash/StateEqual (and StorageHash/StorageEqual) are used for ankerl
// sets/maps.

// ── Dynamic storage type ───────────────────────────────────────────────────
inline bool usesStorage64();

template <typename S> inline S toStorageKey(StateType s);
template <> inline uint64_t toStorageKey<uint64_t>(StateType s) { return s.u.lo; }
template <> inline StateType toStorageKey<StateType>(StateType s) { return s; }

template <typename S> inline S encodeForStorage(StateType state, Transition tr);
template <> inline uint64_t encodeForStorage<uint64_t>(StateType state, Transition tr) {
  return state.u.lo | ((uint64_t)tr.packed << 56);
}
template <> inline StateType encodeForStorage<StateType>(StateType state, Transition tr) {
  return encodeTransition(state, tr);
}

template <typename S> inline StateType storageToState(S s);
template <> inline StateType storageToState<uint64_t>(uint64_t s) { return StateType(s & 0x00FFFFFFFFFFFFFFull, 0); }
template <> inline StateType storageToState<StateType>(StateType s) { return s & STATE_MASK; }

template <typename S> inline Transition storageTrans(S s);
template <> inline Transition storageTrans<uint64_t>(uint64_t s) {
  Transition t;
  t.packed = static_cast<uint8_t>(s >> 56);
  return t;
}
template <> inline Transition storageTrans<StateType>(StateType s) { return decodeTransition(s); }

template <typename S> struct StorageHash;
template <> struct StorageHash<uint64_t> {
  using is_avalanching = void;
  size_t operator()(uint64_t x) const noexcept {
    return size_t(ankerl::unordered_dense::detail::wyhash::hash(x & 0x00FFFFFFFFFFFFFFULL));
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

inline bool usesStorage64() { return g_numPos <= 16; }

inline int countTowersLeft(StateType s) {
  int towersLeft = 0;
  for (int pos = 0; pos < g_numPos; pos++) {
    if (stateGet(s, pos) != 0)
      towersLeft++;
  }
  return towersLeft;
}

// Apply a single permutation (used in non-hot paths and for diagnostics).
inline StateType applyPerm(StateType s, const int *perm) {
  StateType result = 0;
  StateType remaining = s;
  for (int i = 0; i < g_numPos; i++) {
    result = stateSet(result, perm[i], stateGet(remaining, 0));
    remaining = remaining >> 3;
  }
  return result;
}

// Choose the symmetry whose transformed state is lexicographically smallest
// when states are compared by tower type at canonical positions 0..g_numPos-1.
// Always 4 or 8 symmetries; branches unrolled for 4 and 8.
inline int chooseBestSymmetry(StateType s) {
  if (g_numSymm == 4) {
    int cand[4] = {0, 1, 2, 3};
    int candCount = 4;
    for (int canonPos = 0; canonPos < g_numPos && candCount > 1; ++canonPos) {
      int bestVal = 6;
      for (int idx = 0; idx < 4; ++idx) {
        if (idx >= candCount)
          continue;
        int v = stateGet(s, g_symmPermInv[cand[idx]][canonPos]);
        if (v < bestVal)
          bestVal = v;
      }
      int write = 0;
      for (int idx = 0; idx < 4; ++idx) {
        if (idx >= candCount)
          continue;
        int k = cand[idx];
        if (stateGet(s, g_symmPermInv[k][canonPos]) == bestVal)
          cand[write++] = k;
      }
      candCount = write;
    }
    return cand[0];
  }
  // 8 symmetries
  int cand[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  int candCount = 8;
  for (int canonPos = 0; canonPos < g_numPos && candCount > 1; ++canonPos) {
    int bestVal = 6;
    for (int idx = 0; idx < 8; ++idx) {
      if (idx >= candCount)
        continue;
      int v = stateGet(s, g_symmPermInv[cand[idx]][canonPos]);
      if (v < bestVal)
        bestVal = v;
    }
    int write = 0;
    for (int idx = 0; idx < 8; ++idx) {
      if (idx >= candCount)
        continue;
      int k = cand[idx];
      if (stateGet(s, g_symmPermInv[k][canonPos]) == bestVal)
        cand[write++] = k;
    }
    candCount = write;
  }
  return cand[0];
}

inline StateType canonical(StateType s) {
  int bestK = chooseBestSymmetry(s);
  StateType result = StateType(0);
  for (int origPos = 0; origPos < g_numPos; ++origPos) {
    int v = stateGet(s, origPos);
    int canonPos = g_symmPerm[bestK][origPos];
    result = stateSet(result, canonPos, v);
  }
  return result;
}

inline std::pair<StateType, int> canonicalAndChildPos(StateType s, int pos) {
  int bestK = chooseBestSymmetry(s);
  StateType best = StateType(0);
  for (int origPos = 0; origPos < g_numPos; ++origPos) {
    int v = stateGet(s, origPos);
    int canonPos = g_symmPerm[bestK][origPos];
    best = stateSet(best, canonPos, v);
  }
  int canonChildPos = g_symmPerm[bestK][pos];
  return {best, canonChildPos};
}

inline uint32_t neighborPresent(StateType state, int pos) {
  uint32_t present = 0;
  const int *shifts = g_neighborShifts[pos];
  int count = g_neighborCount[pos];
  for (int i = 0; i < count; i++)
    present |= 1u << stateGet(state, shifts[i] / 3);
  return present;
}

inline bool canPlace(uint32_t present, int type) {
  if (type <= 1)
    return true;
  return (present & g_requiredMask[type]) == g_requiredMask[type];
}
