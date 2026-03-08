#pragma once

#include <iostream>

#include "config.h"
#include "grid_config.h"
#include "state_utils.h"

inline const string &bgColor(int t, bool placed) {
  return placed ? COLORS[t] : EMPTY_BG;
}
inline const string &fgColor(int t, bool placed) {
  return placed ? FG_COLORS[t] : EMPTY_FG;
}

// ── Grid display ─────────────────────────────────────────────────────────────

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

