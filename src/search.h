#pragma once

#include "config.h"
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

// High-level solver entrypoint for a single level and unlocked-type
void runSingleLevel();

// When non-null and open, baseline build path is written here instead of cout
extern std::ofstream *g_baselineFile;

// Result of a search run (BFS or sliding window).
struct SearchResult {
  StateType bestState;
  int bestScore;
  long long explored;
  size_t uniqueStates;
};

// Score = sum of position types (0..5). Used by BFS and sliding window.
int scoreGrid(StateType state);

// Print a build path (from reverse BFS). Path is from empty board (all 0) to goal.
void printBuildPathSteps(const std::vector<std::pair<int, int>> &buildPath, std::ostream &out = std::cout);
