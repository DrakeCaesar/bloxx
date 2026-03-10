#!/usr/bin/env python3
"""
Read a best_states file (e.g. level3_type5_best_states.txt) and display each
canonical state as a grid with ANSI colors (same as C++ printGrid).

Usage:
  python scripts/view_best_states.py [file]
  python scripts/view_best_states.py level3_type5_best_states.txt

On Windows, ensure the console uses UTF-8 (e.g. chcp 65001) so half-blocks display.
"""

import re
import sys
from pathlib import Path

# On Windows, force UTF-8 stdout so half-block chars (▐ ▌) and ANSI codes work
if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except AttributeError:
        import io
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", line_buffering=True)

# ANSI colors matching config.h / display.h (index 0 = empty, 1-5 = tower types)
RESET = "\033[0m"
BLACK_TEXT = "\033[30m"
# Background colors for cells
COLORS = [
    "\033[100m",  # Grey - type 0 (empty)
    "\033[44m",   # Blue - type 1
    "\033[41m",   # Red - type 2
    "\033[42m",   # Green - type 3
    "\033[43m",   # Yellow - type 4
    "\033[47m",   # White - type 5
]
# Foreground colors for half-block separators
FG_COLORS = [
    "\033[90m",   # Grey - type 0
    "\033[34m",   # Blue - type 1
    "\033[31m",   # Red - type 2
    "\033[32m",   # Green - type 3
    "\033[33m",   # Yellow - type 4
    "\033[37m",   # White - type 5
]

# Grid patterns from grid_config.h: '.' = valid position, ' ' = hole
LEVEL_PATTERNS = {
    1: [" .. ", "....", "....", " .. "],
    2: ["....", "....", "....", "...."],
    3: [".....", ".....", "  .  ", ".....", "....."],
    4: [
        "......",
        "..  ..",
        " .... ",
        " .... ",
        "..  ..",
        "......",
    ],
    5: [
        "    ...",
        " ......",
        " ......",
        " .. .. ",
        "...... ",
        "...... ",
        "...    ",
    ],
}


def get_positions(level: int) -> tuple[list[tuple[int, int]], int, int]:
    """Return (positions, rows, cols) for the level. positions = [(r,c), ...] in row-major order of '.' cells."""
    pattern = LEVEL_PATTERNS.get(level, LEVEL_PATTERNS[1])
    rows = len(pattern)
    cols = max(len(row) for row in pattern)
    positions = []
    for r, row in enumerate(pattern):
        for c in range(len(row)):
            if row[c] == ".":
                positions.append((r, c))
    return positions, rows, cols


def parse_filename(path: str) -> int | None:
    """Extract level from filename like level3_type5_best_states.txt. Returns None if not matched."""
    name = Path(path).name
    m = re.match(r"level(\d+)_type\d+_best_states\.txt", name, re.IGNORECASE)
    return int(m.group(1)) if m else None


def load_best_states(path: str) -> tuple[int, int, list[str]]:
    """Load file; return (score, count, list of state strings)."""
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    if not lines or not lines[0].startswith("score "):
        raise SystemExit(f"Invalid format: expected 'score N' first line in {path}")
    score = int(lines[0].split()[1])
    if len(lines) < 2 or not lines[1].startswith("count "):
        raise SystemExit(f"Invalid format: expected 'count M' second line in {path}")
    count = int(lines[1].split()[1])
    state_lines = [s for s in lines[2:] if s.strip()]
    if len(state_lines) != count:
        raise SystemExit(
            f"File says count={count} but found {len(state_lines)} state lines in {path}"
        )
    return score, count, state_lines


def state_line_to_grid(state_line: str, positions: list[tuple[int, int]], rows: int, cols: int) -> list[list[int | None]]:
    """Build a rows x cols grid. None = hole, 0-5 = tower type. Position i in state_line maps to positions[i]."""
    grid = [[None] * cols for _ in range(rows)]
    for i, (r, c) in enumerate(positions):
        if i < len(state_line):
            grid[r][c] = int(state_line[i])
    return grid


def exists(grid: list[list[int | None]], r: int, c: int) -> bool:
    """True if (r, c) is a valid cell (not a hole)."""
    if r < 0 or r >= len(grid) or c < 0 or c >= len(grid[0]):
        return False
    return grid[r][c] is not None


def get_type(grid: list[list[int | None]], r: int, c: int) -> int:
    """Cell type 0-5; only call when exists()."""
    t = grid[r][c]
    return t if t is not None else 0


def print_grid(grid: list[list[int | None]], title: str = "") -> None:
    """Print one grid with ANSI colors and half-block separators (same as C++ printGrid)."""
    if title:
        print(title)
    rows, cols = len(grid), len(grid[0])
    for r in range(rows):
        line_parts = []
        # First cell (c=0)
        if exists(grid, r, 0):
            t = get_type(grid, r, 0)
            line_parts.append(f"{FG_COLORS[t]}▐{RESET}")
        else:
            line_parts.append(" ")
        for c in range(cols):
            if exists(grid, r, c):
                t = get_type(grid, r, c)
                line_parts.append(f"{COLORS[t]}{BLACK_TEXT}{t}{RESET}")
            else:
                line_parts.append(" ")
            # Separator between c and c+1
            if c < cols - 1:
                curr_exists = exists(grid, r, c)
                next_exists = exists(grid, r, c + 1)
                if curr_exists and next_exists:
                    t = get_type(grid, r, c)
                    next_t = get_type(grid, r, c + 1)
                    line_parts.append(f"{FG_COLORS[next_t]}{COLORS[t]}▐{RESET}")
                elif curr_exists and not next_exists:
                    t = get_type(grid, r, c)
                    line_parts.append(f"{FG_COLORS[t]}▌{RESET}")
                elif not curr_exists and next_exists:
                    next_t = get_type(grid, r, c + 1)
                    line_parts.append(f"{FG_COLORS[next_t]}▐{RESET}")
                else:
                    line_parts.append(" ")
        # End of row
        if exists(grid, r, cols - 1):
            last_t = get_type(grid, r, cols - 1)
            line_parts.append(f"{FG_COLORS[last_t]}▌{RESET}")
        else:
            line_parts.append(" ")
        print("".join(line_parts))


def main() -> None:
    path = sys.argv[1] if len(sys.argv) > 1 else None
    if not path or path in ("-h", "--help"):
        print(__doc__.strip())
        print("\nUsage: python scripts/view_best_states.py <best_states_file>")
        sys.exit(0 if path in ("-h", "--help") else 1)

    path = Path(path)
    if not path.is_file():
        raise SystemExit(f"File not found: {path}")

    level = parse_filename(str(path))
    if level is None:
        level = 3  # default guess; state line length will be checked
    positions, rows, cols = get_positions(level)

    score, count, state_lines = load_best_states(str(path))

    if not state_lines:
        print(f"Score {score}, no states in file.")
        return

    npos = len(positions)
    print(f"Level (from filename): {level}")
    print(f"Score: {score}  Count: {count}  Positions: {npos}")
    print()

    for idx, line in enumerate(state_lines):
        if len(line) != npos:
            raise SystemExit(
                f"State line {idx + 1}: expected {npos} digits, got {len(line)}. "
                f"Check that the file matches level {level}."
            )
        grid = state_line_to_grid(line, positions, rows, cols)
        print_grid(grid, title=f"--- State {idx + 1}/{count} ---")
        print()


if __name__ == "__main__":
    main()
