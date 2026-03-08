#pragma once

#include "config.h"
#include <fstream>
#include <streambuf>
#include <utility>
#include <vector>

// High-level solver entrypoint for a single level and unlocked-type
void runSingleLevel();

// Baseline output redirection hooks used by the solver
extern std::ofstream *g_baselineFile;
extern std::streambuf *g_originalCoutBuf;

// Result of a search run (BFS, DFS, or sliding window).
struct SearchResult {
  StateType bestState;
  int bestScore;
  long long explored;
  size_t uniqueStates;
};

// Print a build path (from reverse BFS or DFS). Uses PATH_FROM_EMPTY / expand.
void printBuildPathSteps(const std::vector<std::pair<int, int>> &buildPath);

