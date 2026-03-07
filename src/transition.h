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
// Memory: 2 bytes vs. 32 bytes (16-byte-aligned StateType + padding) — ~16x.
struct Transition {
  uint16_t packed; // bits 0-5: child_pos (0-62, 63=sentinel), bits 6-8: prev_type

  Transition(int8_t child_pos = -1, int8_t prev_type = 0)
      : packed((child_pos == -1 ? 63 : child_pos) | (prev_type << 6)) {}

  int8_t getPos() const { return (packed & 0x3F) == 63 ? -1 : (packed & 0x3F); }
  int8_t getType() const { return (packed >> 6) & 0x7; }
};

