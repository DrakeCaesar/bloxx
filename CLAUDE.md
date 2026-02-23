# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Tower Blox Deluxe Solver - A C++ program that finds optimal solutions for the Tower Blox puzzle game. The solver uses state-space search with various algorithms to find the highest-scoring tower configurations.

## Build Commands

```bash
# Configure build (from repository root)
cmake -B build

# Build the solver
cmake --build build --config Release

# Run the solver
./build/tower_solver

# Clean build
rm -rf build && cmake -B build && cmake --build build --config Release
```

The project uses CMake with C++17 standard. Release builds are optimized with `-O3 -march=native` (or MSVC equivalents).

## Core Architecture

### State Representation
The solver uses bit-packing to compress grid states into compact integer types:
- **64-bit states** (uint64_t): For grids with ≤21 positions (levels 1-3)
- **128-bit states** (__uint128_t): For larger grids (levels 4-5)
- Each position stores a tower type (0-4) in 3 bits
- States are canonicalized using grid symmetries to reduce search space

### Key Components

**grid.h** - Grid configuration and utilities:
- `GridConfig`: Stores grid geometry, neighbor relationships, and position constraints
- `createGridFromPattern()`: Creates grids from ASCII patterns
- `createLevel1Grid()` through `createLevel5Grid()`: Pre-defined level configurations
- `Transition`: Compact struct storing move information (parent state, position, new type)
- `STATE_TYPE`: Selected at compile-time based on `LEVEL` constant

**main.cpp** - Search algorithms and state management:
- `buildLookupTables()`: Precomputes neighbor bit-shifts for performance
- `buildSymmetries()`: Auto-detects valid grid symmetries (8 dihedral transforms)
- `canonical()`: Returns minimum state across all symmetry transformations
- `neighborPresent()` / `canPlace()`: Core game logic for placement validation

### Search Algorithms

The solver supports multiple search modes (selected via `SEARCH_MODE` constant):

1. **BFS (0)**: Breadth-first search, finds optimal solutions but memory-intensive
2. **DFS (1)**: Depth-first search, lower memory but may not find optimal solutions
3. **IDDFS (2)**: Iterative deepening DFS, balances memory and optimality
4. **Sliding Window (3)**: Novel optimization for large grids, processes 3x3 windows in zigzag order

### Sliding Window Optimization

The sliding window mode (`SEARCH_MODE = 3`) is designed for large grids:
- Processes grid in overlapping kernel windows (configurable patterns)
- Uses zigzag traversal from 4 corners (A→B→C→D passes)
- Reconstructs build path via reverse deconstruction + BFS scaffolding
- Handles invalid steps by finding intermediate "scaffolding" moves

## Game Rules

Tower placement follows dependency rules:
- Type 0 (Blue): Can be placed anywhere (no requirements)
- Type 1 (Red): Requires at least one Type 0 neighbor
- Type 2 (Green): Requires at least one Type 0 AND one Type 1 neighbor
- Type 3 (Yellow): Requires Types 0, 1, and 2 neighbors
- Type 4 (White): Requires all lower types (0-3) as neighbors

Score is sum of all tower types on the grid (max type = 4 per position).

## Configuration

Key constants in grid.h:
- `LEVEL`: 1-5, selects grid configuration
- `SEARCH_MODE`: 0-3, algorithm selection
- `MAX_UNLOCKED_TYPE`: 0-4, maximum tower type available
- `MAX_DEPTH`: BFS depth limit per sliding window
- `SCAFFOLD_DEPTH`: BFS depth for path reconstruction scaffolding
- `KERNEL_PATTERNS`: Window shapes for sliding window mode
- `KERNEL_INDEX`: Selects which kernel pattern to use

## Performance Considerations

- **Symmetry reduction**: Auto-detects and applies valid dihedral symmetries (D4 group)
- **State canonicalization**: All states stored in canonical form to avoid duplicates
- **Lookup tables**: Precomputed neighbor bit-shifts eliminate runtime calculations
- **Robin-hood hashing**: Used for visited state tracking (fast, cache-friendly)
- **Memory usage**: Grows with visited states; BFS mode uses most memory

## Development Notes

When modifying the solver:
- The `LEVEL` constant must be a `#define` for preprocessor selection of state type
- State packing assumes max 4 tower types (fits in 3 bits) - changing this requires bit-packing logic updates
- Path reconstruction is only native for BFS mode; other modes use post-hoc BFS
- The sliding window optimization trades absolute optimality for scalability on large grids
- Window kernel patterns affect convergence - different shapes work better for different grid geometries
