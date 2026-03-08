#include "search_dfs.h"
#include "search_impl.h"
#include "search.h"
#include "display.h"
#include "config.h"
#include "transition.h"

#include <robin_hood.h>
#include <climits>
#include <iostream>
#include <vector>

using namespace std;

// Convert state sequence (path from start to goal) to build steps (pos, newType).
static vector<pair<int, int>>
stateSequenceToSteps(const vector<StateType> &path) {
  vector<pair<int, int>> steps;
  steps.reserve(path.empty() ? 0 : (path.size() - 1));
  for (size_t i = 0; i + 1 < path.size(); i++) {
    StateType s = path[i], t = path[i + 1];
    int stepTypeMin = PATH_FROM_EMPTY ? 0 : 1;
    bool found = false;
    for (int pos = 0; pos < g_numPos && !found; pos++) {
      for (int newType = stepTypeMin; newType <= 5; newType++) {
        if (stateGet(s, pos) == newType)
          continue;
        if (canonical(stateSet(s, pos, newType)) == t) {
          steps.push_back({pos, newType});
          found = true;
          break;
        }
      }
    }
    if (!found)
      steps.push_back({-1, 0});
  }
  return steps;
}

// Path = states only; parent = path[depth-1]. We only skip a state if we've
// already reached it at a strictly shorter depth. When we improve (better
// score or same score with shorter path) we store the full path.
static void dfsVisit(
    StateType cur, int curScore, int depth, vector<StateType> &path,
    robin_hood::unordered_map<StateType, int, StateHash, StateEqual> &minDepth,
    SearchResult &result, vector<StateType> &bestPath, int &bestPathLength) {
  result.explored++;
  int pathLen = (int)path.size() - 1;
  if (curScore > result.bestScore ||
      (curScore == result.bestScore && pathLen < bestPathLength)) {
    result.bestScore = curScore;
    result.bestState = cur;
    bestPathLength = pathLen;
    bestPath = path;
    cout << "New best score: " << result.bestScore << " (path length "
         << bestPathLength << ")\n";
    printGrid(cur);
  }
  if (depth >= DFS_MAX_DEPTH)
    return;
  for (int pos = 0; pos < g_numPos; pos++) {
    int curType = stateGet(cur, pos);
    uint32_t present = neighborPresent(cur, pos);
    for (int newType = 0; newType <= g_maxType[pos]; newType++) {
      if (!FORWARD_TRANSITIONS[curType][newType])
        continue;
      if (!canPlace(present, newType))
        continue;
      auto [next, childPos] =
          canonicalAndChildPos(stateSet(cur, pos, newType), pos);
      int nextDepth = depth + 1;
      auto it = minDepth.find(next);
      if (it != minDepth.end() && it->second <= nextDepth)
        continue;
      minDepth[next] = nextDepth;
      path.push_back(next);
      dfsVisit(next, curScore + (newType - curType), nextDepth, path, minDepth,
               result, bestPath, bestPathLength);
      path.pop_back();
    }
  }
}

void runDFSSearch(SearchResult &result) {
  robin_hood::unordered_map<StateType, int, StateHash, StateEqual> minDepth;
  minDepth.reserve(1ull << 22);
  vector<StateType> path;
  vector<StateType> bestPath;
  int bestPathLength = INT_MAX;

  StateType startState = PATH_FROM_EMPTY ? StateType(0) : g_baseState;
  minDepth[startState] = 0;
  path.push_back(startState);

  result.bestState = startState;
  result.bestScore = 0;
  result.explored = 0;

  dfsVisit(startState, 0, 0, path, minDepth, result, bestPath, bestPathLength);
  result.uniqueStates = minDepth.size();

  if (g_baselineFile && g_baselineFile->is_open()) {
    g_originalCoutBuf = cout.rdbuf();
    cout.rdbuf(g_baselineFile->rdbuf());
  }
  vector<pair<int, int>> buildPath = stateSequenceToSteps(bestPath);
  if (!buildPath.empty())
    printBuildPathSteps(buildPath);
  if (g_originalCoutBuf) {
    cout.rdbuf(g_originalCoutBuf);
    g_originalCoutBuf = nullptr;
  }
}
