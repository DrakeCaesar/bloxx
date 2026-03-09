#pragma once

#include <cstdint>

// ── Path reconstruction transition ───────────────────────────────────────────
// Stores only (child_pos, prev_type) — no parent state pointer needed.
//
// child_pos  = position index in the *child's* canonical state that changed.
// prev_type  = the type at child_pos in the parent (what we revert to).
//
// Parent recovery: canonical(stateSet(child, child_pos, prev_type))
// Proof: stateSet(child, child_pos, prev_type) == applyPerm(parent, k) for the
// symmetry k that mapped parent→child, so canonicalizing gives back parent. ✓
//
// Encoding: packed = child_pos * 6 + prev_type  (0–203)
//   34 positions × 6 types = 204 distinct values ≤ 255 → fits in one byte.
//   255 = sentinel (start/final state, no transition).
//
// Memory: 1 byte vs. 32 bytes (16-byte-aligned StateType + padding) — ~32x.
struct Transition {
  uint8_t packed; // child_pos*6 + prev_type (0-203), 255 = sentinel

  Transition(int8_t child_pos = -1, int8_t prev_type = 0)
      : packed(child_pos == -1 ? 255 : static_cast<uint8_t>(child_pos * 6 + prev_type)) {}

  int8_t getPos() const { return packed == 255 ? -1 : packed / 6; }
  int8_t getType() const { return packed == 255 ? 0 : packed % 6; }
};
