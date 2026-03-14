#!/usr/bin/env python3
"""
Read all best_states files (best_states/levelN_typeU.txt), count type transitions
from paths (state_0 -> state_1 -> ... -> state_n), and print 6x6 tables plus
config.h-style matrices. Types: 0=Empty, 1=Blue, 2=Red, 3=Green, 4=Yellow, 5=White.

Usage:
  python scripts/count_transitions.py [--all-paths] [--dir DIR]
  --all-paths: use every path (default: only shortest paths per file)
  --dir: directory containing best_states files (default: best_states)
"""

import argparse
import re
from pathlib import Path

TYPE_NAMES = ["Empty", "Blue", "Red", "Green", "Yellow", "White"]
CONFIG_NAMES = ["empty", "blue", "red", "green", "yellow", "white"]


def load_best_states(path: Path) -> tuple[int, int, list[tuple[str, list[str] | None]]]:
    """Load file; return (score, count, list of (state_str, path_sequence)).
    path_sequence is None if no path, else [state_0, ..., state_n] (each npos digits)."""
    text = path.read_text()
    lines = [line.rstrip("\n") for line in text.splitlines()]
    if not lines or not lines[0].startswith("score "):
        raise SystemExit(f"Invalid format: expected 'score N' first line in {path}")
    score = int(lines[0].split()[1])
    if len(lines) < 2 or not lines[1].startswith("count "):
        raise SystemExit(f"Invalid format: expected 'count M' second line in {path}")
    count = int(lines[1].split()[1])
    raw_lines = [s for s in lines[2:] if s.strip()]
    if len(raw_lines) != count:
        raise SystemExit(
            f"File says count={count} but found {len(raw_lines)} state lines in {path}"
        )
    entries: list[tuple[str, list[str] | None]] = []
    for line in raw_lines:
        parts = line.split()
        if not parts:
            entries.append(("", None))
            continue
        state_str = parts[0]
        path_list: list[str] | None = None
        if len(parts) >= 2:
            try:
                path_len = int(parts[1])
                need = 2 + path_len + 1  # path_len+1 states
                if path_len >= 0 and len(parts) >= need:
                    path_list = parts[2:need]
            except ValueError:
                pass
        entries.append((state_str, path_list))
    return score, count, entries


def transitions_from_path(path_list: list[str]) -> list[tuple[int, int]]:
    """From path [state_0, state_1, ..., state_n], return list of (from_type, to_type)
    for each step (the single cell that changed)."""
    out: list[tuple[int, int]] = []
    n = len(path_list[0])
    for i in range(len(path_list) - 1):
        a, b = path_list[i], path_list[i + 1]
        if len(a) != n or len(b) != n:
            continue
        for pos in range(n):
            if a[pos] != b[pos]:
                from_t = int(a[pos])
                to_t = int(b[pos])
                if 0 <= from_t <= 5 and 0 <= to_t <= 5:
                    out.append((from_t, to_t))
                break
    return out


def find_best_states_files(directory: Path) -> list[Path]:
    """Return sorted list of best_states files (levelN_typeU.txt or *_best_states.txt)."""
    pattern = re.compile(r"level(\d+)_type(\d+)\.txt$", re.IGNORECASE)
    pattern2 = re.compile(r"level(\d+)_type(\d+)_best_states\.txt$", re.IGNORECASE)
    found: list[tuple[int, int, Path]] = []
    for p in directory.iterdir():
        if not p.is_file() or p.suffix.lower() != ".txt":
            continue
        name = p.name
        m = pattern.match(name) or pattern2.match(name)
        if m:
            found.append((int(m.group(1)), int(m.group(2)), p))
    found.sort(key=lambda t: (t[0], t[1]))
    return [p for _, _, p in found]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Count type transitions from best_states path files and print 6x6 tables."
    )
    parser.add_argument(
        "--all-paths",
        action="store_true",
        help="Use all paths; default is only shortest paths per file",
    )
    repo_root = Path(__file__).resolve().parent.parent
    parser.add_argument(
        "--dir",
        type=Path,
        default=repo_root / "best_states",
        metavar="DIR",
        help="Directory containing levelN_typeU.txt files (default: repo best_states/)",
    )
    args = parser.parse_args()
    use_shortest_only = not args.all_paths

    best_states_dir = args.dir
    if not best_states_dir.is_dir():
        raise SystemExit(f"Directory not found: {best_states_dir}")

    files = find_best_states_files(best_states_dir)
    if not files:
        raise SystemExit(f"No best_states files found in {best_states_dir}")

    # Global 6x6 count: cnt[from][to]
    cnt: list[list[int]] = [[0] * 6 for _ in range(6)]

    for path in files:
        score, count, entries = load_best_states(path)
        paths_with_len: list[tuple[list[str], int]] = []
        for _state_str, path_list in entries:
            if path_list is None or len(path_list) < 2:
                continue
            paths_with_len.append((path_list, len(path_list) - 1))  # steps = len-1
        if not paths_with_len:
            continue
        min_steps = min(steps for _, steps in paths_with_len)
        for path_list, steps in paths_with_len:
            if use_shortest_only and steps != min_steps:
                continue
            for from_t, to_t in transitions_from_path(path_list):
                cnt[from_t][to_t] += 1

    # Print forward table
    print("Forward (from -> to)")
    print()
    print("%-8s" % "from\\to", end="")
    for j in range(6):
        print("%8s" % TYPE_NAMES[j], end="")
    print()
    print("-" * 56)
    for i in range(6):
        print("%-8s" % TYPE_NAMES[i], end="")
        for j in range(6):
            print("%8s" % (cnt[i][j] if cnt[i][j] > 0 else "."), end="")
        print()
    print("-" * 56)
    print("%-8s" % "TOTAL", end="")
    for j in range(6):
        col = sum(cnt[i][j] for i in range(6))
        print("%8s" % (col if col > 0 else "."), end="")
    print()

    # Backward table (cell [i][j] = forward [j][i])
    print()
    print("Backward (from -> to)")
    print()
    print("%-8s" % "from\\to", end="")
    for j in range(6):
        print("%8s" % TYPE_NAMES[j], end="")
    print()
    print("-" * 56)
    for i in range(6):
        print("%-8s" % TYPE_NAMES[i], end="")
        for j in range(6):
            val = cnt[j][i]
            print("%8s" % (val if val > 0 else "."), end="")
        print()
    print("-" * 56)
    print("%-8s" % "TOTAL", end="")
    for j in range(6):
        col = sum(cnt[j][i] for i in range(6))
        print("%8s" % (col if col > 0 else "."), end="")
    print()

    # config.h paste blocks
    config_names_padded = [n.ljust(7) for n in CONFIG_NAMES]
    print()
    print("--- Paste into config.h: FORWARD_TRANSITIONS ---")
    print()
    for i in range(6):
        row = (
            "    /* "
            + config_names_padded[i]
            + "*/ {"
            + ", ".join("1" if cnt[i][j] > 0 else "0" for j in range(6))
            + "},"
        )
        print(row)
    print()
    print("--- Paste into config.h: BACKWARD_TRANSITIONS ---")
    print()
    for i in range(6):
        row = (
            "    /* "
            + config_names_padded[i]
            + "*/ {"
            + ", ".join("1" if cnt[j][i] > 0 else "0" for j in range(6))
            + "},"
        )
        print(row)


if __name__ == "__main__":
    main()
