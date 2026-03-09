#include "search_bfs.h"
#include "search_impl.h"
#include "search.h"
#include "config.h"
#include "display.h"

#include <ankerl/unordered_dense.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <climits>
#include <deque>
#include <iostream>
#include <mutex>
#include <optional>
#include <thread>
#include <vector>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif

using QueueEntry = std::pair<StateType, int>;

#ifdef _WIN32
static void trimWorkingSetIfHighMemory() {
  MEMORYSTATUSEX mem = {};
  mem.dwLength = sizeof(mem);
  if (!GlobalMemoryStatusEx(&mem))
    return;
  if (mem.dwMemoryLoad > 80)
    SetProcessWorkingSetSize(GetCurrentProcess(), (SIZE_T)-1, (SIZE_T)-1);
}
#endif

template <typename S> struct ReverseExpansion {
  S next;
  int nextTowers;
};

template <typename ParentSet>
static void pruneSideImpl(ParentSet &parentMap, std::deque<QueueEntry> &queue,
                          bool keepLowest) {
  using S = typename ParentSet::value_type;
  if (parentMap.size() <= REVERSE_BFS_PRUNE_THRESHOLD || queue.empty())
    return;

  int minTowers = g_numPos;
  int maxTowers = 0;
  for (const auto &entry : queue) {
    int tc = entry.second;
    minTowers = std::min(minTowers, tc);
    maxTowers = std::max(maxTowers, tc);
  }

  int cutoff = keepLowest ? (minTowers + REVERSE_BFS_PRUNE_LOW_DELTA)
                          : (maxTowers - REVERSE_BFS_PRUNE_HIGH_DELTA);
  std::vector<QueueEntry> toProcess;
  toProcess.reserve(queue.size());
  for (const auto &entry : queue) {
    if ((keepLowest && entry.second <= cutoff) ||
        (!keepLowest && entry.second >= cutoff)) {
      toProcess.push_back(entry);
    }
  }

  if (toProcess.empty()) {
    if (keepLowest)
      toProcess.push_back(
          *std::min_element(queue.begin(), queue.end(),
                            [](const QueueEntry &a, const QueueEntry &b) {
                              return a.second < b.second;
                            }));
    else
      toProcess.push_back(
          *std::max_element(queue.begin(), queue.end(),
                            [](const QueueEntry &a, const QueueEntry &b) {
                              return a.second > b.second;
                            }));
  }

  ParentSet pruned;
  pruned.reserve(toProcess.size() * 4);
  for (const auto &entry : toProcess)
    pruned.emplace(*parentMap.find(toStorageKey<S>(entry.first)));

  std::deque<QueueEntry> worklist(toProcess.begin(), toProcess.end());
  while (!worklist.empty()) {
    StateType cur = worklist.front().first;
    worklist.pop_front();
    auto it = parentMap.find(toStorageKey<S>(cur));
    if (it == parentMap.end())
      continue;
    auto tr = storageTrans<S>(*it);
    if (tr.getPos() == -1)
      continue;
    StateType parent = canonical(stateSet(cur, tr.getPos(), tr.getType()));
    if (parentMap.count(toStorageKey<S>(parent)) &&
        pruned.find(toStorageKey<S>(parent)) == pruned.end()) {
      pruned.emplace(*parentMap.find(toStorageKey<S>(parent)));
      worklist.push_back({parent, countTowersLeft(parent)});
    }
  }

  parentMap.swap(pruned);
  queue.assign(toProcess.begin(), toProcess.end());
}

template <typename S>
static void
reverseExpandChunkImpl(const QueueEntry *start, const QueueEntry *end,
                       std::vector<ReverseExpansion<S>> *out,
                       StateType pathStartState, int reverseBaseType) {
  const auto &backTable =
      PATH_FROM_EMPTY ? BACKWARD_TRANSITIONS_EMPTY : BACKWARD_TRANSITIONS;
  for (; start != end; ++start) {
    StateType cur = start->first;
    int curTowers = start->second;
    if (cur == pathStartState)
      continue;
    for (int pos = 0; pos < g_numPos; pos++) {
      int curType = stateGet(cur, pos);
      uint32_t present = neighborPresent(cur, pos);
      if (!canPlace(present, curType))
        continue;
      for (int prevType = 0; prevType <= g_maxType[pos]; prevType++) {
        if (!backTable[curType][prevType])
          continue;
        auto [next, childPos] =
            canonicalAndChildPos(stateSet(cur, pos, prevType), pos);
        int nextTowers = curTowers - (curType != reverseBaseType ? 1 : 0) +
                         (prevType != reverseBaseType ? 1 : 0);
        out->push_back({encodeForStorage<S>(next, Transition{(int8_t)childPos,
                                                             (int8_t)curType}),
                        nextTowers});
      }
    }
  }
}

template <typename S> struct ReverseTaskContext {
  std::deque<QueueEntry> *qReverse;
  ankerl::unordered_dense::segmented_set<S, StorageHash<S>, StorageEqual<S>>
      *parentReverse;
  StateType pathStartState;
  std::mutex *coutMutex;
  std::atomic<bool> *reverseFinished;
  size_t *totalReverseVisited;
  int reverseBaseType;
};

template <typename S> static void runReverseTask(ReverseTaskContext<S> *ctx) {
  std::deque<QueueEntry> &qReverse = *ctx->qReverse;
  auto &parentReverse = *ctx->parentReverse;
  StateType pathStartState = ctx->pathStartState;
  std::mutex &coutMutex = *ctx->coutMutex;
  std::atomic<bool> &reverseFinished = *ctx->reverseFinished;
  size_t &totalReverseVisited = *ctx->totalReverseVisited;
  int reverseBaseType = ctx->reverseBaseType;

  size_t depth = 0;
  while (!qReverse.empty()) {
    size_t layerCount = qReverse.size();
    {
      std::lock_guard<std::mutex> lock(coutMutex);
      std::cout << "[Reverse] Layer " << depth << ": expanding "
                << formatNum(layerCount) << " states ("
                << std::max(1, std::min(REVERSE_BFS_NUM_THREADS, (int)layerCount))
                << " threads)..." << std::endl;
    }
    auto layerStart = std::chrono::steady_clock::now();

    std::vector<QueueEntry> layer;
    layer.reserve(layerCount);
    for (size_t i = 0; i < layerCount; i++) {
      layer.push_back(qReverse.front());
      qReverse.pop_front();
    }

    int minTowersLeft = g_numPos;
    int maxTowersLeft = 0;
    bool reachedBase = false;
    for (const auto &entry : layer) {
      if (entry.first == pathStartState) {
        reachedBase = true;
        break;
      }
      minTowersLeft = std::min(minTowersLeft, entry.second);
      maxTowersLeft = std::max(maxTowersLeft, entry.second);
    }
    if (reachedBase) {
      qReverse.clear();
      auto layerEnd = std::chrono::steady_clock::now();
      double expandMs =
          std::chrono::duration_cast<std::chrono::milliseconds>(layerEnd -
                                                                layerStart)
              .count();
      std::lock_guard<std::mutex> lock(coutMutex);
      std::cout << "[Reverse] Layer " << depth << " done in "
                << formatTimeMs((long long)expandMs) << " (expand "
                << formatTimeMs((long long)expandMs) << "). Reached base state."
                << std::endl;
      break;
    }

    const int numThreads =
        std::max(1, std::min(REVERSE_BFS_NUM_THREADS, (int)layer.size()));
    const size_t chunkSize =
        (layer.size() + (size_t)numThreads - 1) / (size_t)numThreads;
    std::vector<std::vector<ReverseExpansion<S>>> results((size_t)numThreads);
    for (int c = 0; c < numThreads; c++)
      results[c].reserve(layer.size() * 2);

    if (numThreads == 1) {
      reverseExpandChunkImpl<S>(layer.data(), layer.data() + layer.size(),
                                &results[0], pathStartState, reverseBaseType);
    } else {
      std::vector<std::thread> revThreads;
      for (int c = 0; c < numThreads; c++) {
        const QueueEntry *chunkStart = layer.data() + (size_t)c * chunkSize;
        const QueueEntry *chunkEnd =
            layer.data() +
            std::min((size_t)(c + 1) * chunkSize, layer.size());
        if (chunkStart < chunkEnd)
          revThreads.emplace_back(reverseExpandChunkImpl<S>, chunkStart,
                                  chunkEnd, &results[c], pathStartState,
                                  reverseBaseType);
      }
      for (auto &t : revThreads)
        t.join();
    }

    {
      size_t inserted = 0;
      constexpr size_t trimCheckInterval = 50000;
      for (int c = 0; c < numThreads; c++) {
        for (const auto &r : results[c]) {
          if (parentReverse.emplace(r.next).second) {
            totalReverseVisited++;
            qReverse.push_back({storageToState<S>(r.next), r.nextTowers});
            if (++inserted % trimCheckInterval == 0) {
#ifdef _WIN32
              trimWorkingSetIfHighMemory();
#endif
            }
          }
        }
      }
    }

    auto layerEnd = std::chrono::steady_clock::now();
    double expandMs =
        std::chrono::duration_cast<std::chrono::milliseconds>(layerEnd -
                                                              layerStart)
            .count();
    double pruneMs = 0.0;
    if (!qReverse.empty()) {
      auto pruneStart = std::chrono::steady_clock::now();
      pruneSideImpl(parentReverse, qReverse, true);
      pruneMs = std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::steady_clock::now() - pruneStart)
                    .count();
    }
    double totalMs = expandMs + pruneMs;
    {
      std::lock_guard<std::mutex> lock(coutMutex);
      std::cout << "[Reverse] Layer " << depth << " done in "
                << formatTimeMs((long long)totalMs) << " (expand "
                << formatTimeMs((long long)expandMs) << ", prune "
                << formatTimeMs((long long)pruneMs) << "). Queue size now "
                << formatNum(qReverse.size()) << ", map size "
                << formatNum(parentReverse.size()) << "." << std::endl;
    }
#ifdef _WIN32
    trimWorkingSetIfHighMemory();
#endif
    depth++;
  }
  {
    std::lock_guard<std::mutex> lock(coutMutex);
    std::cout << "[Reverse] Queue empty. Reverse BFS finished." << std::endl;
  }
  reverseFinished = true;
}

template <typename S> struct ForwardTaskContext {
  std::deque<QueueEntry> *qForward;
  ankerl::unordered_dense::segmented_set<S, StorageHash<S>, StorageEqual<S>>
      *parentForward;
  StateType pathStartState;
  int reverseBaseType;
  std::atomic<bool> *reverseFinished;
  std::mutex *coutMutex;
};

template <typename S> static void runForwardTask(ForwardTaskContext<S> *ctx) {
  std::deque<QueueEntry> &qForward = *ctx->qForward;
  auto &parentForward = *ctx->parentForward;
  int reverseBaseType = ctx->reverseBaseType;
  std::atomic<bool> &reverseFinished = *ctx->reverseFinished;
  std::mutex &coutMutex = *ctx->coutMutex;

  size_t forwardDepth = 0;
  while (!qForward.empty()) {
    if (reverseFinished) {
      std::lock_guard<std::mutex> lock(coutMutex);
      std::cout << "[Forward] Reverse finished; stopping immediately."
                << std::endl;
      break;
    }
    size_t forwardLayerCount = qForward.size();
    {
      std::lock_guard<std::mutex> lock(coutMutex);
      std::cout << "[Forward] Layer " << forwardDepth << ": expanding "
                << formatNum(forwardLayerCount) << " states..." << std::endl;
    }
    auto forwardLayerStart = std::chrono::steady_clock::now();
    int minForwardTowersLeft = g_numPos;
    int maxForwardTowersLeft = 0;
    bool anyForwardProcessed = false;
    bool stoppedByReverse = false;

    for (size_t i = 0; i < forwardLayerCount; i++) {
      if (reverseFinished) {
        stoppedByReverse = true;
        break;
      }
      StateType cur = qForward.front().first;
      int curTowers = qForward.front().second;
      qForward.pop_front();
      if (!anyForwardProcessed) {
        minForwardTowersLeft = curTowers;
        maxForwardTowersLeft = curTowers;
        anyForwardProcessed = true;
      } else {
        minForwardTowersLeft = std::min(minForwardTowersLeft, curTowers);
        maxForwardTowersLeft = std::max(maxForwardTowersLeft, curTowers);
      }

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
          if (parentForward
                  .emplace(encodeForStorage<S>(
                      next, Transition{(int8_t)childPos, (int8_t)curType}))
                  .second) {
            int nextTowers = curTowers - (curType != reverseBaseType ? 1 : 0) +
                             (newType != reverseBaseType ? 1 : 0);
            qForward.push_back({storageToState<S>(next), nextTowers});
          }
        }
      }
    }

    auto forwardLayerEnd = std::chrono::steady_clock::now();
    double expandMs = std::chrono::duration_cast<std::chrono::milliseconds>(
                          forwardLayerEnd - forwardLayerStart)
                          .count();
    double pruneMs = 0.0;
    if (!stoppedByReverse && !qForward.empty()) {
      auto pruneStart = std::chrono::steady_clock::now();
      pruneSideImpl(parentForward, qForward, false);
      pruneMs = std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::steady_clock::now() - pruneStart)
                    .count();
    }
    double totalMs = expandMs + pruneMs;
    {
      std::lock_guard<std::mutex> lock(coutMutex);
      std::cout << "[Forward] Layer " << forwardDepth << " done in "
                << formatTimeMs((long long)totalMs) << " (expand "
                << formatTimeMs((long long)expandMs) << ", prune "
                << formatTimeMs((long long)pruneMs) << "). Queue size now "
                << formatNum(qForward.size()) << ", map size "
                << formatNum(parentForward.size());
      if (stoppedByReverse)
        std::cout << " (stopped mid-layer: reverse finished)";
      std::cout << "." << std::endl;
    }
    if (stoppedByReverse)
      break;
    forwardDepth++;
  }
  {
    std::lock_guard<std::mutex> lock(coutMutex);
    std::cout << "[Forward] Forward BFS finished (queue empty or reverse "
                 "finished)."
              << std::endl;
  }
}

template <typename S>
static int forwardDepthOfImpl(
    const ankerl::unordered_dense::segmented_set<S, StorageHash<S>, StorageEqual<S>>
        &parentForward,
    StateType s) {
  int d = 0;
  StateType cur = s;
  while (parentForward.count(toStorageKey<S>(cur))) {
    auto tr = storageTrans<S>(*parentForward.find(toStorageKey<S>(cur)));
    if (tr.getPos() == -1)
      break;
    cur = canonical(stateSet(cur, tr.getPos(), tr.getType()));
    d++;
  }
  return d;
}

template <typename S>
static int reverseDepthOfImpl(
    const ankerl::unordered_dense::segmented_set<S, StorageHash<S>, StorageEqual<S>>
        &parentReverse,
    StateType s) {
  int d = 0;
  StateType cur = s;
  while (parentReverse.count(toStorageKey<S>(cur))) {
    auto tr = storageTrans<S>(*parentReverse.find(toStorageKey<S>(cur)));
    if (tr.getPos() == -1)
      break;
    cur = canonical(stateSet(cur, tr.getPos(), tr.getType()));
    d++;
  }
  return d;
}

template <typename S>
static std::vector<std::pair<int, int>>
buildPathFromReverseBfsImpl(StateType finalState) {
  const StateType pathStartState = PATH_FROM_EMPTY ? StateType(0) : g_baseState;
  const int reverseBaseType = PATH_FROM_EMPTY ? 0 : 1;

#ifdef _WIN32
  SetProcessWorkingSetSize(GetCurrentProcess(), (SIZE_T)-1, (SIZE_T)-1);
#endif

  ankerl::unordered_dense::segmented_set<S, StorageHash<S>, StorageEqual<S>>
      parentReverse;
  ankerl::unordered_dense::segmented_set<S, StorageHash<S>, StorageEqual<S>>
      parentForward;

  std::deque<QueueEntry> qReverse;
  std::deque<QueueEntry> qForward;
  if (RUN_PARALLEL_FORWARD_BFS) {
    parentForward.emplace(encodeForStorage<S>(pathStartState, {}));
    qForward.push_back({pathStartState, countTowersLeft(pathStartState)});
  }
  parentReverse.emplace(encodeForStorage<S>(finalState, {}));
  qReverse.push_back({finalState, countTowersLeft(finalState)});

  std::mutex coutMutex;
  std::atomic<bool> reverseFinished{false};
  size_t totalReverseVisited = 0;

  {
    std::lock_guard<std::mutex> lock(coutMutex);
    if (RUN_PARALLEL_FORWARD_BFS)
      std::cout << "[Main] Starting reverse BFS (to completion) and forward "
                   "BFS in parallel. Forward will stop as soon as reverse "
                   "finishes."
                << std::endl;
    else
      std::cout << "[Main] Reverse BFS only (parallel forward disabled)."
                << std::endl;
  }

  ReverseTaskContext<S> revCtx{
      &qReverse,        &parentReverse,       pathStartState, &coutMutex,
      &reverseFinished, &totalReverseVisited, reverseBaseType};
  std::thread tReverse(runReverseTask<S>, &revCtx);
  std::thread tForward;
  std::optional<ForwardTaskContext<S>> fwdCtx;
  if (RUN_PARALLEL_FORWARD_BFS) {
    fwdCtx =
        ForwardTaskContext<S>{&qForward,       &parentForward,   pathStartState,
                             reverseBaseType, &reverseFinished, &coutMutex};
    tForward = std::thread(runForwardTask<S>, &*fwdCtx);
  }
  tReverse.join();
  if (RUN_PARALLEL_FORWARD_BFS) {
    {
      std::lock_guard<std::mutex> lock(coutMutex);
      std::cout << "[Main] Reverse finished. Not waiting for forward to "
                   "complete; stopping forward and finding best intercept with "
                   "current state."
                << std::endl;
    }
    tForward.join();
  }

  std::cout << "\nReverse BFS visited " << formatNum(totalReverseVisited)
            << " states (map size after pruning: "
            << formatNum(parentReverse.size()) << ").\n";

  if (g_baselineFile && g_baselineFile->is_open()) {
    g_originalCoutBuf = std::cout.rdbuf();
    std::cout.rdbuf(g_baselineFile->rdbuf());
  }

  std::vector<std::pair<int, int>> steps;

  if (!RUN_PARALLEL_FORWARD_BFS) {
    if (parentReverse.count(toStorageKey<S>(pathStartState)) == 0) {
      if (g_originalCoutBuf) {
        std::cout.rdbuf(g_originalCoutBuf);
        g_originalCoutBuf = nullptr;
      }
      std::lock_guard<std::mutex> lock(coutMutex);
      std::cout << "[Main] Reverse BFS did not reach path start. Cannot build "
                   "path."
                << std::endl;
      return {};
    }
    StateType cur = pathStartState;
    while (cur != finalState) {
      auto tr = storageTrans<S>(*parentReverse.find(toStorageKey<S>(cur)));
      StateType parent = canonical(stateSet(cur, tr.getPos(), tr.getType()));
      bool foundStep = false;
      int stepTypeMin = PATH_FROM_EMPTY ? 0 : 1;
      for (int p = 0; p < g_numPos && !foundStep; p++) {
        for (int t = stepTypeMin; t <= g_maxType[p] && !foundStep; t++) {
          if (t == stateGet(cur, p))
            continue;
          if (canonical(stateSet(cur, p, t)) == parent) {
            steps.push_back({p, t});
            foundStep = true;
          }
        }
      }
      cur = parent;
    }
    return steps;
  }

  StateType meetState = pathStartState;
  int bestTotal = INT_MAX;
  size_t interceptCount = 0;
  for (const auto &elem : parentForward) {
    StateType s = storageToState<S>(elem);
    if (parentReverse.count(toStorageKey<S>(s)) == 0)
      continue;
    interceptCount++;
    int fd = forwardDepthOfImpl<S>(parentForward, s);
    int rd = reverseDepthOfImpl<S>(parentReverse, s);
    int total = fd + rd;
    if (total < bestTotal) {
      bestTotal = total;
      meetState = s;
    }
  }

  {
    std::lock_guard<std::mutex> lock(coutMutex);
    std::cout << "[Main] Found " << formatNum(interceptCount)
              << " states in both maps. Best intercept has forward depth "
              << forwardDepthOfImpl<S>(parentForward, meetState)
              << ", reverse depth "
              << reverseDepthOfImpl<S>(parentReverse, meetState) << ", total "
              << formatNum((long long)bestTotal) << " steps." << std::endl;
  }

  if (interceptCount == 0) {
    if (g_originalCoutBuf) {
      std::cout.rdbuf(g_originalCoutBuf);
      g_originalCoutBuf = nullptr;
    }
    std::lock_guard<std::mutex> lock(coutMutex);
    std::cout << "[Main] No intercept: no state reached by both searches. "
                 "Cannot reconstruct path."
              << std::endl;
    return {};
  }

  std::vector<std::pair<int, int>> stepsForward;
  {
    StateType cur = meetState;
    while (true) {
      auto tr = storageTrans<S>(*parentForward.find(toStorageKey<S>(cur)));
      if (tr.getPos() == -1)
        break;
      int childPos = tr.getPos();
      int newType = stateGet(cur, childPos);
      StateType parent = canonical(stateSet(cur, childPos, tr.getType()));
      int posInParent = -1;
      for (int p = 0; p < g_numPos; p++) {
        if (stateGet(parent, p) == newType)
          continue;
        if (canonical(stateSet(parent, p, newType)) == cur) {
          posInParent = p;
          break;
        }
      }
      if (posInParent < 0)
        posInParent = childPos;
      stepsForward.push_back({posInParent, newType});
      cur = parent;
    }
  }
  std::reverse(stepsForward.begin(), stepsForward.end());

  {
    StateType cur = meetState;
    while (true) {
      auto tr = storageTrans<S>(*parentReverse.find(toStorageKey<S>(cur)));
      if (tr.getPos() == -1)
        break;
      StateType parent = canonical(stateSet(cur, tr.getPos(), tr.getType()));
      bool foundStep = false;
      int stepTypeMin = PATH_FROM_EMPTY ? 0 : 1;
      for (int p = 0; p < g_numPos && !foundStep; p++) {
        for (int t = stepTypeMin; t <= g_maxType[p] && !foundStep; t++) {
          if (t == stateGet(cur, p))
            continue;
          if (canonical(stateSet(cur, p, t)) == parent) {
            steps.push_back({p, t});
            foundStep = true;
          }
        }
      }
      cur = parent;
    }
  }

  steps.insert(steps.begin(), stepsForward.begin(), stepsForward.end());
  return steps;
}

std::vector<std::pair<int, int>> buildPathFromReverseBfs(StateType finalState) {
  if (usesStorage64())
    return buildPathFromReverseBfsImpl<uint64_t>(finalState);
  else
    return buildPathFromReverseBfsImpl<StateType>(finalState);
}
