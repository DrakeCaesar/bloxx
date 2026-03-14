# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Tower Blox Deluxe Solver - A C++20 program that finds optimal solutions for the Tower Blox puzzle game through state-space search. The solver uses bit-packing, symmetry reduction, and novel sliding window optimizations to efficiently search large state spaces.

## Build Commands

```bash
# Configure build (from repository root)
cmake -B build

# Build the solver
cmake --build build --config Release

# Run the solver (Windows)
./build/Release/tower_solver.exe
# or
./build/tower_solver  (MinGW/Linux)

# Clean build
rm -rf build && cmake -B build && cmake --build build --config Release

# Build for specific target (Visual Studio)
cmake --build build --config Release --target tower_solver
```

The project uses CMake with C++20 standard. Single dependency: [ankerl::unordered_dense](https://github.com/martinus/unordered_dense) (FetchContent). Release builds optimized with `/O3 /GL /arch:AVX2` (MSVC) or `-O3 -march=native -flto` (GCC/Clang).

## Core Architecture

### File Structure

**driver.cpp** - Entry point and baseline generation:
- Orchestrates multi-level/multi-type baseline generation
- Creates output files in format `level{N}_type{U}_baseline.txt`
- Configures console UTF-8 mode on Windows

**search.cpp** - BFS/sliding-window orchestration, symmetry, lookup tables:
- State management, symmetry handling, sliding window BFS
- Hot-path functions marked `inline` for performance
- Global lookup tables (`g_neighborShifts`, `g_symmPerm`) for zero-overhead access

**search_bfs.cpp** - Reverse BFS for path reconstruction (mode 0 and after sliding window).

**Headers** (provide modularity and separation of concerns):
- **config.h** - State type, constants, compile-time configuration (level/unlocked type ranges)
- **grid_config.h** - Grid creation, neighbor computation, max type analysis; uses `ankerl::unordered_dense::map` for coord→index
- **state_utils.h** - Bit manipulation (`stateSet`, `stateGet`)
- **transition.h** - Compact 16-bit transition encoding for path reconstruction
- **search.h** - Solver entry point, global file redirection
- **display.h** - Unicode grid rendering with ANSI colors

### State Representation

- **SIMD-based 128-bit state**: Custom `uint128_simd` type using `__m128i` union
  - Stable ABI on Windows (MSVC/clang), passed in XMM registers
  - Supports all bitwise ops, shifts, and comparisons
  - Handles up to 42 positions × 3 bits/position=126 bits
- **State type**: 128-bit for all levels (multi-level support). Single-level could use 64-bit for levels 1–3 via `getCurrentLevel()`.

### Performance-Critical Optimizations

**Lookup Tables** (precomputed in `buildLookupTables()`):
- `g_neighborShifts[49][4]` - Pre-multiplied neighbor indices × 3 bits
- `g_neighborCount[49]` - Neighbor count per position
- `g_maxType[49]` - Maximum achievable tower type per position
- Padded to fixed size (49×4) to eliminate vector pointer indirection

**Symmetry Reduction** (auto-detected in `buildSymmetries()`):
- Tests all 8 dihedral transforms (identity, 90°, 180°, 270°, mirror variants)
- Only keeps symmetries where every position maps to another valid position
- Identity always first (index 0) for early exit
- Canonical form = minimum state across all valid symmetries

**State Canonicalization**:
- All states stored and compared in canonical form
- Eliminates duplicate search paths due to grid symmetries
- Reduces search space by 2-8x depending on grid

### Search Algorithms

Three search modes (configure in config.h via `SEARCH_MODE`):

1. **BFS (0)**: Breadth-first search using queue
 - Guaranteed optimal solutions
 - Memory-intensive (stores all visited states)
 - Path reconstruction via reverse BFS in search_bfs.cpp

2. **Sliding Window (1)**: Novel optimization for large grids
 - Default and most effective for levels 3-5
 - See detailed section below
 - Path reconstruction via reverse BFS in search_bfs.cpp

3. **Sliding window only (2)**: Same sliding window as mode 1, but no path reconstruction
 - Runs sliding window to find optimal board state only; no reverse BFS to empty/start

### Sliding Window Optimization (SEARCH_MODE=1)

**Core Algorithm**:
- Processes grid in overlapping windows following a JPEG-like zigzag path
- 4 corner-based passes (A→B→C→D) cover entire grid
- Each window: BFS with depth limit `MAX_DEPTH` (default 10)
- Window kernel configurable via `KERNEL_PATTERNS` (cross, square, diamond, etc.)

**Zigzag Traversal**:
- Builds diagonal-wise path within each window
- Alternates direction per diagonal (true zigzag)
- Example 3×3: (0,0)→(0,1)→(1,0)→(2,0)→(1,1)→(0,2)→(1,2)→(2,1)→(2,2)

**Path Reconstruction**:
1. Reverse deconstruction from final state by removing towers in reverse order
2. Validate each step is legal (precomputed during forward search)
3. If invalid step found:
   - Run BFS from valid state with limit `SCAFFOLD_DEPTH` (default 10)
   - Find intermediate "scaffolding" moves to bridge the gap
   - Continue with corrected path

**Threading**: Reverse BFS uses `REVERSE_BFS_NUM_THREADS` (default 8) threads per layer for parallel expansion

### Placement Rules & Game Logic

Tower types are numbered 1-5 (0 = empty/unplaced). Path finding starts from empty (all 0).

Tower placement follows dependency rules verified by `canPlace()`:
- **Type 1 (Blue)**: Can be placed anywhere (no requirements)
- **Type 2 (Red)**: Requires at least one Type 1 neighbor
- **Type 3 (Green)**: Requires at least one Type 1 AND one Type 2 neighbor
- **Type 4 (Yellow)**: Requires Types 1, 2, and 3 neighbors
- **Type 5 (White)**: Requires all lower types (1-4) as neighbors

Score = sum of position types for all positions (type 0 = 0, type 1 = 1, …, type 5 = 5)

`neighborPresent()` computes bitmask once per (state, position) for efficiency

## Configuration

**In config.h** - Key compile-time constants:

**Grid & Search**:
- `SEARCH_MODE`: Algorithm selector (0=BFS, 1=Sliding Window, 2=Sliding window only)
- `KERNEL_INDEX`: Selects window kernel pattern (0-5, see patterns below)

**Sliding Window**:
- `MAX_DEPTH`: BFS depth limit per window (default: 10)
- `SCAFFOLD_DEPTH`: BFS depth for scaffolding moves during path reconstruction (default: 10)

**Multi-Level Baseline** (single source of truth for level and tower type):
- `LEVEL_START`, `LEVEL_END`: Level range to process (inclusive), e.g. 1–5
- `UNLOCKED_TYPE_START`, `UNLOCKED_TYPE_END`: Tower type range to test (e.g. 1–5)
- `REVERSE_BFS_NUM_THREADS`: Threads for parallel reverse BFS (default: 8)

**Window Kernel Patterns** (`KERNEL_PATTERNS` array):
```cpp
0: Cross shape       (5 positions)    " . ", "...", " . "
1: 3×3 square        (9 positions)    "...", "...", "..."
2: Diamond/Plus      (13 positions)   "  .  ", " ... ", ".....", " ... ", "  .  "
3: 5×5 square        (25 positions)
4: Plus shape        (9 positions)
5: L-shape           (5 positions)
```

## Multi-Level Baseline Generation

**Workflow**:
- Driver processes `LEVEL_START` to `LEVEL_END` (e.g., levels 1-5)
- For each level, tests all unlocked tower types from `UNLOCKED_TYPE_START` to `UNLOCKED_TYPE_END`
- Generates 25 baseline files: `level{N}_type{U}_baseline.txt`

**Baseline File Format**:
- Console output from reverse BFS section
- Valid build sequence for verification
- Used for regression testing and solution verification

**Running Single-Level**:
- Set `LEVEL_START = LEVEL_END = desired_level`
- Set `UNLOCKED_TYPE_START = UNLOCKED_TYPE_END = desired_type`
- Run solver: baseline saved, console shows sliding window progress

## Performance Characteristics

**Scaling**:
- **Levels 1-2**: BFS finds optimal solutions quickly (small state space)
- **Level 3**: Sliding window recommended (moderate complexity)
- **Levels 4-5**: Sliding window essential (large state space)

**Memory Usage**:
- Large visited/parent containers use `ankerl::unordered_dense::segmented_set` / `segmented_map` for lower peak memory and smoother growth.
- BFS: Grows with visited states (GBs for large grids)
- Sliding window: Window-sized frontier (MBs to low GBs)

**CPU**:
- Lookup tables eliminate runtime computation in hot path
- SIMD 128-bit ops efficiently vectorized by compiler
- Threading only during reverse BFS (path reconstruction)

**Hashing**:
- All hash containers use [ankerl::unordered_dense](https://github.com/martinus/unordered_dense). Custom hashes (StateHash, StorageHash, `hash<std::pair<int,int>>`) use wyhash and `is_avalanching` for high-quality hashing without extra mixing.

## Development Guidelines

**When Modifying**:
- **State type**: 128-bit (`StateType = StateType128`) for multi-level support. For single-level, level is from `getCurrentLevel()` (set by driver).
- **Bit packing**: 3 bits/position × 42 positions = 126 bits max (values 0-5; 0=empty, 1-5=tower types). Changing tower types requires updating all shift/mask logic.
- **Grid patterns**: Use `createGridFromPattern()` with '.' for valid positions
- **New levels**: Add to switch statement in `createGrid()` and recompute baselines
- **Symmetry**: Always canonicalize states before storing/looking up
- **Hot path**: `neighborPresent()`, `canPlace()`, and inner loops must remain `inline`

**Testing**:
- Generate baselines for new configurations
- Compare against previous baselines for regression detection
- Verify path reconstruction produces valid sequences

**Common Tasks**:
- **Optimize kernel**: Adjust `KERNEL_INDEX` and regenerate baselines for comparison
- **Tune depth**: Modify `MAX_DEPTH` and `SCAFFOLD_DEPTH` (tradeoff: speed vs. solution quality)
- **Profile**: Use optimized Release build; hotspots in state canonicalization and neighbor checks
