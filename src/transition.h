#pragma once

#include "config.h"

// ── Path reconstruction transition ───────────────────────────────────────────
// Transition struct optimized for board size
// For grids ≤52 positions: pack pos (6 bits) + type (3 bits) + sentinel into
// uint16_t. This reduces Transition size and improves cache utilization.
struct Transition {
  StateType parent;
  uint16_t packed; // bits 0-5: pos (0-63, 63=sentinel), bits 6-8: type (0-5)

  Transition(StateType p = 0, int8_t pos = -1, int8_t type = 0)
      : parent(p), packed((pos == -1 ? 63 : pos) | (type << 6)) {}

  int8_t getPos() const { return (packed & 0x3F) == 63 ? -1 : (packed & 0x3F); }
  int8_t getType() const { return (packed >> 6) & 0x7; }
};

