#pragma once

#include <iostream>
#include <ostream>
#include <algorithm>
#include <vector>
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

inline std::string padLeft(const std::string &s, size_t width, char fill = ' ') {
  if (s.size() >= width)
    return s;
  return std::string(width - s.size(), fill) + s;
}
inline std::string padRight(const std::string &s, size_t width, char fill = ' ') {
  if (s.size() >= width)
    return s;
  return s + std::string(width - s.size(), fill);
}

// Fixed width for reverse BFS status labels so progress bars align in one column.
inline constexpr size_t REVERSE_BFS_LABEL_WIDTH = 42;

inline size_t utf8Codepoints(const std::string &s) {
  // Counts UTF-8 codepoints (assumes valid UTF-8). Good enough for console
  // column width when using single-cell glyphs (incl. Private Use Area).
  size_t n = 0;
  for (size_t i = 0; i < s.size(); ) {
    unsigned char c = (unsigned char)s[i];
    if ((c & 0x80) == 0x00) i += 1;
    else if ((c & 0xE0) == 0xC0) i += 2;
    else if ((c & 0xF0) == 0xE0) i += 3;
    else if ((c & 0xF8) == 0xF0) i += 4;
    else i += 1;
    n++;
  }
  return n;
}

// Progress bar: "[=========>    ] 45% (9/20)" when total > 0; width in chars
inline std::string progressBar(size_t current, size_t total, int width = 30) {
  if (total == 0) return "[  no total  ] 0% (0/0)";
  if (current > total) current = total;
  size_t pct = total ? (current * 100) / total : 0;

  // Private Use Area glyphs (as shown in the screenshot): EE00-EE05.
  // Empty:  EE00 (left cap), EE01 (mid), EE02 (right cap)
  // Filled: EE03 (left cap), EE04 (mid), EE05 (right cap)
  const std::string L_EMPTY = "\xEE\xB8\x80";     // U+EE00
  const std::string MID_EMPTY = "\xEE\xB8\x81";   // U+EE01
  const std::string R_EMPTY = "\xEE\xB8\x82";     // U+EE02
  const std::string L_FILL = "\xEE\xB8\x83";      // U+EE03
  const std::string MID_FILL = "\xEE\xB8\x84";    // U+EE04
  const std::string R_FILL = "\xEE\xB8\x85";      // U+EE05

  int inner = (std::max)(1, width);
  int filled = total ? (int)((current * (size_t)inner) / total) : 0;
  if (filled < 0) filled = 0;
  if (filled > inner) filled = inner;
  std::string bar;
  bar.reserve((size_t)inner * 3 + 64);
  bar += (filled > 0 ? L_FILL : L_EMPTY);
  for (int i = 0; i < inner; i++)
    bar += (i < filled ? MID_FILL : MID_EMPTY);
  bar += (filled >= inner ? R_FILL : R_EMPTY);
  bar += " ";
  bar += padLeft(std::to_string((unsigned)pct), 3);
  bar += "% (";
  std::string totalS = formatNum(total);
  std::string curS = formatNum(current);
  bar += padLeft(curS, totalS.size());
  bar += "/";
  bar += totalS;
  bar += ")";
  return bar;
}

// Phase header for clear "what's running" (single line)
inline void printPhaseHeader(const std::string &title) {
  std::cout << "\n═══ " << title << " ═══\n";
}
inline void printOperation(const std::string &operation, const std::string &detail) {
  std::cout << "  - " << operation;
  if (!detail.empty()) std::cout << "  " << detail;
  std::cout << "\n";
}
inline void printPhaseTiming(const std::string &phaseName, long long elapsedMs) {
  std::cout << "  [time] " << phaseName << ": " << formatTimeMs(elapsedMs) << "\n";
}

// ── Table helpers (build borders/rows from column widths) ───────────────────

enum class TableAlign { Left, Right, Center };

inline std::string repeatUtf8(const std::string &unit, int count) {
  if (count <= 0)
    return {};
  std::string out;
  out.reserve((size_t)count * unit.size());
  for (int i = 0; i < count; i++)
    out += unit;
  return out;
}

inline std::string fitCell(const std::string &s, int width, TableAlign align) {
  if (width <= 0)
    return {};
  if ((int)s.size() >= width)
    return s.substr(0, (size_t)width);
  int pad = width - (int)s.size();
  if (align == TableAlign::Left)
    return s + std::string((size_t)pad, ' ');
  if (align == TableAlign::Right)
    return std::string((size_t)pad, ' ') + s;
  int left = pad / 2;
  int right = pad - left;
  return std::string((size_t)left, ' ') + s + std::string((size_t)right, ' ');
}

struct TableStyle {
  std::string h = "─";
  std::string v = "│";
  std::string tl = "┌", tm = "┬", tr = "┐";
  std::string ml = "├", mm = "┼", mr = "┤";
  std::string bl = "└", bm = "┴", br = "┘";
};

inline void printTableBorderTop(const std::vector<int> &widths,
                                const TableStyle &st = TableStyle{}) {
  std::cout << st.tl;
  for (size_t i = 0; i < widths.size(); i++) {
    std::cout << repeatUtf8(st.h, widths[i]);
    std::cout << (i + 1 == widths.size() ? st.tr : st.tm);
  }
  std::cout << "\n";
}

inline void printTableBorderMid(const std::vector<int> &widths,
                                const TableStyle &st = TableStyle{}) {
  std::cout << st.ml;
  for (size_t i = 0; i < widths.size(); i++) {
    std::cout << repeatUtf8(st.h, widths[i]);
    std::cout << (i + 1 == widths.size() ? st.mr : st.mm);
  }
  std::cout << "\n";
}

inline void printTableBorderBottom(const std::vector<int> &widths,
                                   const TableStyle &st = TableStyle{}) {
  std::cout << st.bl;
  for (size_t i = 0; i < widths.size(); i++) {
    std::cout << repeatUtf8(st.h, widths[i]);
    std::cout << (i + 1 == widths.size() ? st.br : st.bm);
  }
  std::cout << "\n";
}

inline void printTableRow(const std::vector<std::string> &cells,
                          const std::vector<int> &widths,
                          const std::vector<TableAlign> &aligns,
                          const TableStyle &st = TableStyle{}) {
  std::cout << st.v;
  for (size_t i = 0; i < widths.size(); i++) {
    const std::string cell = i < cells.size() ? cells[i] : std::string{};
    const TableAlign a = i < aligns.size() ? aligns[i] : TableAlign::Left;
    std::cout << fitCell(cell, widths[i], a) << st.v;
  }
  std::cout << "\n";
}

// Live-updating single-line status (redraws in-place using '\r').
// IMPORTANT: Only use when the output stream is a terminal (not redirected).
struct DynamicStatusLine {
  bool enabled = false;
  size_t lastCols = 0;
  std::string lastLine;

  explicit DynamicStatusLine(bool enable) : enabled(enable) {}

  void clear() {
    if (!enabled)
      return;
    std::cout << '\r';
    if (lastCols)
      std::cout << std::string(lastCols, ' ');
    std::cout << '\r' << std::flush;
    lastCols = 0;
    lastLine.clear();
  }

  void update(const std::string &line) {
    if (!enabled) {
      std::cout << line << "\n";
      return;
    }
    const size_t cols = utf8Codepoints(line);
    std::cout << '\r' << line;
    if (cols < lastCols)
      std::cout << std::string(lastCols - cols, ' ');
    std::cout << std::flush;
    lastCols = (std::max)(lastCols, cols);
    lastLine = line;
  }

  void redraw() {
    if (!enabled || lastLine.empty())
      return;
    update(lastLine);
  }

  void finish(const std::string &finalLine = "") {
    if (!enabled) {
      if (!finalLine.empty())
        std::cout << finalLine << "\n";
      return;
    }
    if (!finalLine.empty())
      update(finalLine);
    std::cout << "\n";
    lastCols = 0;
    lastLine.clear();
  }
};

inline DynamicStatusLine *g_activeStatusLine = nullptr;

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

