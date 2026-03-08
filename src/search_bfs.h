#pragma once

#include "config.h"
#include <utility>
#include <vector>

// Reverse BFS from final state back to path start. Used by sliding window for
// path reconstruction.
std::vector<std::pair<int, int>> buildPathFromReverseBfs(StateType finalState);
