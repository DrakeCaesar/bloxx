#pragma once

#include "config.h"

// ════════════════════════════════════════════════════════════════════════════
// State packing utilities
// ════════════════════════════════════════════════════════════════════════════
// Each position stores a value 0-5 (0=empty, 1-5=tower types, fits in 3 bits).
// Up to 21 positions × 3 bits = 63 bits → fits in uint64_t.

inline StateType stateSet(StateType s, int pos, int val) {
  s &= ~(StateType(7) << (pos * 3));
  s |= StateType(val) << (pos * 3);
  return s;
}

inline int stateGet(StateType s, int pos) { return (s >> (pos * 3)) & 7; }

