#pragma once

#include "search.h"

// Depth-limited DFS with min-depth map and in-memory path. Used when
// SEARCH_MODE == 2.
void runDFSSearch(SearchResult &result);
