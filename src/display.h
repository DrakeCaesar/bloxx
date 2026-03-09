#pragma once

#include <iostream>
#include <string>

#include "config.h"
#include "grid_config.h"
#include "state_utils.h"

// Format integer with space thousands separators (e.g. 1234567 -> "1 234 567")
inline std::string formatNum(size_t n) {
  std::string s = std::to_string(n);
  for (int i = (int)s.size() - 3; i > 0; i -= 3)
    s.insert(i, 1, ' ');
  return s;
}
inline std::string formatNum(long long n) {
  if (n < 0) return "-" + formatNum((size_t)(-n));
  return formatNum((size_t)n);
}

// Format milliseconds as "M min S s" or "S s" (e.g. 125000 -> "2 min 5 s", 45000 -> "45 s", 500 -> "0.5 s")
inline std::string formatTimeMs(long long ms) {
  if (ms < 0) return "-" + formatTimeMs(-ms);
  if (ms < 1000) {
    // Sub-second: show as "0.5 s" etc.
    std::string s = std::to_string(static_cast<double>(ms) / 1000.0);
    size_t dot = s.find('.');
    if (dot != std::string::npos) {
      while (s.size() > dot + 1 && s.back() == '0') s.pop_back();
      if (s.back() == '.') s.pop_back();
    }
    return s + " s";
  }
  long long sec = ms / 1000;
  long long min = sec / 60;
  sec %= 60;
  if (min > 0)
    return std::to_string(min) + " min " + std::to_string(sec) + " s";
  return std::to_string(sec) + " s";
}

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

