#include <ankerl/unordered_dense.h>
#include <atomic>
#include <chrono>
#include <climits>
#include <deque>
#include <functional>
#include <iostream>
#include <mutex>
#include <robin_hood.h>
#include <thread>
#include <vector>

#include "config.h"
#include "display.h"
#include "grid_config.h"
#include "search.h"
#include "state_utils.h"
#include "transition.h"

using namespace std;

// Global file stream for baseline output (only used for reverse BFS section)
ofstream *g_baselineFile = nullptr;
streambuf *g_originalCoutBuf = nullptr;

// Hash specialization for uint128_simd (must be after including config.h)
namespace robin_hood {
template <> struct hash<uint128_simd> {
  size_t operator()(const uint128_simd &x) const noexcept {
    uint64_t lo = (uint64_t)x;
    uint64_t hi = (uint64_t)(x >> 64);
    return size_t(lo ^ hi);
  }
};
} // namespace robin_hood

// ── Hot-path lookup tables ──────────────────────────────────────────────────
// Flat fixed-size arrays eliminate vector pointer indirection in the inner
// loop. Max grid 7×7 = 49 positions, max 4 orthogonal neighbors each. Neighbor
// bit-shifts are precomputed (neighborIndex * 3) to avoid multiplies in
// canPlace. All positions are padded to exactly 4 entries; unused entries use
// shift 63 which reads bits 63-65 of the uint64_t — always 0 (type 0/empty),
// so they never set any bit that placement checks require (bits 1+).

static int g_numPos;
static int g_neighborShifts[49][4]; // precomputed: neighborIndex * 3
static int g_neighborCount[49];
static int g_maxType[49];
static StateType g_baseState; // all positions set to type 1 (Blue)
// Type T requires neighbors of types 1..T-1: mask = bits 1..T-1 set
static const uint32_t g_requiredMask[6] = {0, 0, 2, 6, 14, 30};

// ── Symmetry ────────────────────────────────────────────────────────────────
// Auto-detects valid symmetries by trying all 8 dihedral transforms against
// actual grid positions. Keeps only those where every position maps to another
// valid position. Identity is always first (index 0).

static int g_numSymm;
static int g_symmPerm[8][49]; // g_symmPerm[s][i] = image of position i under
                              // symmetry s

inline StateType applyPerm(StateType s, const int *perm) {
  StateType result = 0;
  for (int i = 0; i < g_numPos; i++)
    result |= StateType((s >> (i * 3)) & 7) << (perm[i] * 3);
  return result;
}

// Canonical form = minimum StateType value over all symmetry images.
// The minimum is a unique representative of the equivalence class.
inline StateType canonical(StateType s) {
  if (g_numSymm == 1)
    return s; // No symmetries, skip computation
  StateType best = s;
  for (int k = 1; k < g_numSymm; k++) {
    StateType t = applyPerm(s, g_symmPerm[k]);
    if (t < best) {
      best = t;
      if (best == g_baseState)
        return g_baseState; // Early exit: can't get smaller than all-blue
    }
  }
  return best;
}

static void buildSymmetries() {
  int R = GRID.rows, C = GRID.cols;
  // All 8 elements of the dihedral group D4:
  //   identity, 180° rot, 90° CW, 90° CCW, H-flip, V-flip, diagonal,
  //   anti-diagonal
  using T = function<pair<int, int>(int, int)>;
  T transforms[8] = {
      [](int r, int c) { return make_pair(r, c); },
      [R, C](int r, int c) { return make_pair(R - 1 - r, C - 1 - c); },
      [R](int r, int c) { return make_pair(c, R - 1 - r); },
      [C](int r, int c) { return make_pair(C - 1 - c, r); },
      [C](int r, int c) { return make_pair(r, C - 1 - c); },
      [R](int r, int c) { return make_pair(R - 1 - r, c); },
      [](int r, int c) { return make_pair(c, r); },
      [R, C](int r, int c) { return make_pair(C - 1 - c, R - 1 - r); },
  };

  g_numSymm = 0;
  for (int t = 0; t < 8; t++) {
    int perm[49];
    bool valid = true;
    for (int i = 0; i < g_numPos && valid; i++) {
      auto [r, c] = GRID.positions[i];
      auto [nr, nc] = transforms[t](r, c);
      auto it = GRID.coordToIdx.find({nr, nc});
      if (it == GRID.coordToIdx.end())
        valid = false;
      else
        perm[i] = it->second;
    }
    if (valid) {
      for (int i = 0; i < g_numPos; i++)
        g_symmPerm[g_numSymm][i] = perm[i];
      g_numSymm++;
    }
  }
  cout << "Symmetries detected: " << g_numSymm << "  (search space reduced ~"
       << g_numSymm << "x)\n";
}

static void buildLookupTables() {
  g_numPos = GRID.numPositions();
  for (int i = 0; i < g_numPos; i++) {
    const auto &nb = GRID.neighbors[i];
    int count = (int)nb.size();
    g_neighborCount[i] = count;
    g_maxType[i] = GRID.maxType[i];
    for (int j = 0; j < count; j++)
      g_neighborShifts[i][j] = nb[j] * 3;
    for (int j = count; j < 4; j++)
      g_neighborShifts[i][j] = 63; // safe padding: reads always-zero bits
  }
  g_baseState = StateType(0);
  for (int i = 0; i < g_numPos; i++)
    g_baseState = stateSet(g_baseState, i, 1);
}

// ── Placement rule ───────────────────────────────────────────────────────────
// Tower of type T at pos requires one neighbor of each type 1..T-1.
// neighborPresent() computes the bitmask once per (state, pos); canPlace()
// then checks it cheaply for each candidate type.

inline uint32_t neighborPresent(StateType state, int pos) {
  uint32_t present = 0;
  const int *shifts = g_neighborShifts[pos];
  int count = g_neighborCount[pos];
  for (int i = 0; i < count; i++)
    present |= 1u << static_cast<int>((state >> shifts[i]) & 7);
  return present;
}

inline bool canPlace(uint32_t present, int type) {
  if (type <= 1)
    return true;
  return (present & g_requiredMask[type]) == g_requiredMask[type];
}

// ── Path reconstruction helpers ──────────────────────────────────────────────
static int scoreGrid(StateType state) {
  int score = 0;
  for (int i = 0; i < g_numPos; i++)
    score += stateGet(state, i) - 1;
  return score;
}

// expandPathForEmptyBoard: converts a path from the solver's "all-1 = blue"
// model to a path that starts from a truly empty board.
//
// For each step (pos, newType), before we can upgrade pos, we need:
//   1. pos itself to be explicitly placed (as blue) if not already
//   2. For each required neighbor type 1..newType-1, at least one neighbor
//      must be explicitly placed with that type.
//      - Neighbors with type > 1 were placed by earlier explicit steps.
//      - Neighbors with type == 1 are "implicit blue" and need an explicit
//        "place blue" step inserted here.
//
// Tracks a physical (non-canonical) state so neighbor lookups against
// GRID.neighbors[] remain consistent with position indices.
static vector<pair<int, int>>
expandPathForEmptyBoard(const vector<pair<int, int>> &steps) {
  vector<pair<int, int>> out;
  out.reserve(steps.size() * 2);
  uint64_t placed = 0; // bit i = position i has been explicitly placed
  StateType phys = g_baseState; // canonical state — mirrors what printExpandedPath tracks

  for (auto [pos, newType] : steps) {
    // For each required neighbor type req in 1..newType-1,
    // ensure at least one placed neighbor has that type.
    // Unplaced neighbors with type 1 are "implicit blue" in the solver model;
    // insert an explicit place-blue step for one such neighbor if needed.
    for (int req = 1; req < newType; req++) {
      bool satisfied = false;
      int candidateNb = -1;
      for (int nb : GRID.neighbors[pos]) {
        int nbType = stateGet(phys, nb);
        if (nbType == req) {
          if ((placed >> nb) & 1) {
            satisfied = true; // already explicitly placed with correct type
            break;
          } else if (req == 1) {
            // Implicit blue (unplaced, type 1) — can satisfy req=1 once placed
            candidateNb = nb;
          }
        }
      }
      if (!satisfied && candidateNb != -1) {
        out.push_back({candidateNb, 1});
        placed |= (1ull << candidateNb);
        // Blue placements don't change canonical state (1 is already 1)
      }
    }

    // Execute the actual step — keep phys canonical to mirror BFS state space
    out.push_back({pos, newType});
    placed |= (1ull << pos);
    phys = canonical(stateSet(phys, pos, newType));
  }

  return out;
}

static void printExpandedPath(const vector<pair<int, int>> &expandedSteps) {
  StateType canonState = g_baseState;   // canonical state — for validity checks
  StateType displayState = g_baseState; // physical state — stable orientation for display
  uint64_t placedMask = 0;              // physical positions that have been placed

  auto advanceState = [&](int pos, int newType) {
    int usedSymm = 0;
    for (int k = 0; k < g_numSymm; k++) {
      if (applyPerm(displayState, g_symmPerm[k]) == canonState) {
        usedSymm = k;
        break;
      }
    }
    int inversePerm[49];
    for (int i = 0; i < g_numPos; i++)
      inversePerm[g_symmPerm[usedSymm][i]] = i;

    int physPos = inversePerm[pos];
    canonState = canonical(stateSet(canonState, pos, newType));
    displayState = stateSet(displayState, physPos, newType);
    placedMask |= (1ull << physPos);
  };

  // Pre-pass to compute total steps including trailing implicit blues.
  for (auto [pos, newType] : expandedSteps)
    advanceState(pos, newType);

  int missingCount = 0;
  for (int i = 0; i < g_numPos; i++) {
    if (!((placedMask >> i) & 1))
      missingCount++;
  }
  int totalSteps = (int)expandedSteps.size() + missingCount;

  // Reset for actual print pass.
  canonState = g_baseState;
  displayState = g_baseState;
  placedMask = 0;

  cout << "\n--- Build path (" << totalSteps << " steps) ---\n\n";
  cout << "Start (empty board):\n";
  printGrid(g_baseState, uint64_t(0));
  int validSteps = 0, invalidSteps = 0;
  int stepNum = 1;

  for (auto [pos, newType] : expandedSteps) {
    // pos is a position index in canonical space.
    // Find which symmetry maps our displayState to its canonical form,
    // then invert it to find the physical position.
    int usedSymm = 0;
    for (int k = 0; k < g_numSymm; k++) {
      if (applyPerm(displayState, g_symmPerm[k]) == canonState) {
        usedSymm = k;
        break;
      }
    }
    // Build inverse permutation: inversePerm[canonical_pos] = physical_pos
    int inversePerm[49];
    for (int i = 0; i < g_numPos; i++)
      inversePerm[g_symmPerm[usedSymm][i]] = i;

    int physPos = inversePerm[pos];
    bool wasPlaced = (placedMask >> physPos) & 1;
    int oldType = wasPlaced ? stateGet(displayState, physPos) : -1;

    // Validity check in canonical space
    uint32_t present = neighborPresent(canonState, pos);
    bool valid = canPlace(present, newType);
    if (valid)
      validSteps++;
    else
      invalidSteps++;

    // Update canonical state
    canonState = canonical(stateSet(canonState, pos, newType));
    // Update physical display state at the physical position
    displayState = stateSet(displayState, physPos, newType);
    placedMask |= (1ull << physPos);

    auto [r, c] = GRID.positions[physPos];
    cout << "\nStep " << stepNum << ": pos " << physPos << " (" << r << "," << c
         << ")  ";
    if (oldType == -1) {
      cout << "empty -> " << COLORS[newType] << BLACK_TEXT
           << COLOR_NAMES[newType] << RESET << " (" << newType << ")";
    } else {
      cout << COLORS[oldType] << BLACK_TEXT << COLOR_NAMES[oldType] << RESET
           << " (" << oldType << ") -> " << COLORS[newType] << BLACK_TEXT
           << COLOR_NAMES[newType] << RESET << " (" << newType << ")";
    }
    cout << (valid ? "  \033[32m\u2713\033[0m"
                   : "  \033[31m\u2717 INVALID\033[0m")
         << "\n";
    printGrid(displayState, placedMask);
    stepNum++;
  }

  // Any physical position never explicitly placed is an implicit blue that
  // stayed blue the entire game. Show a final step for each one.
  for (int i = 0; i < g_numPos; i++) {
    if (!((placedMask >> i) & 1)) {
      placedMask |= (1ull << i);
      // displayState already has 1 (blue) at this position
      auto [r, c] = GRID.positions[i];
      cout << "\nStep " << stepNum << ": pos " << i << " (" << r << "," << c
           << ")  ";
      cout << "empty -> " << COLORS[1] << BLACK_TEXT << COLOR_NAMES[1] << RESET
           << " (1)";
      cout << "  \033[32m\u2713\033[0m"
           << "\n";
      printGrid(displayState, placedMask);
      stepNum++;
      validSteps++;
    }
  }

  int reachedScore = scoreGrid(displayState);
  cout << "\nBuild path: " << totalSteps << " steps, " << validSteps
       << " valid, " << invalidSteps << " invalid";
  if (invalidSteps == 0) {
    cout << " \033[32m\u2014 ALL VALID\033[0m";
  } else {
    cout << " \033[31m\u2014 HAS INVALID STEPS\033[0m";
  }
  cout << "\nReached score: " << reachedScore << "\n";
}

template <typename MapType>
inline void printBuildPath(StateType bestState, const MapType &visited) {
  vector<pair<int, int>> steps; // (pos, newType)
  StateType cur = bestState;
  while (visited.at(cur).getPos() != -1) {
    const auto &tr = visited.at(cur);
    steps.push_back({tr.getPos(), tr.getType()});
    cur = tr.parent;
  }
  reverse(steps.begin(), steps.end());
  printExpandedPath(expandPathForEmptyBoard(steps));
}

struct SearchResult {
  StateType bestState;
  int bestScore;
  long long explored;
  size_t uniqueStates;
};

// BFS from base state (all blue) stopping as soon as target is reached.
// Returns a parent-pointer map containing the shortest path to target.
inline robin_hood::unordered_map<StateType, Transition>
bfsShortestPath(StateType target) {
  robin_hood::unordered_map<StateType, Transition> path;
  path.reserve(1ull << 24); // ~16M elements
  struct Entry {
    StateType state;
    int score;
  };
  vector<Entry> q;
  size_t head = 0;

  path[g_baseState] = {g_baseState, -1, -1};
  q.push_back({g_baseState, 0});

  while (head < q.size()) {
    auto [cur, curScore] = q[head++];
    if (cur == target)
      break;
    for (int pos = 0; pos < g_numPos; pos++) {
      int curType = stateGet(cur, pos);
      if (curType == 5)
        continue;
      uint32_t present = neighborPresent(cur, pos);
      for (int newType = 1; newType <= g_maxType[pos]; newType++) {
        if (newType == curType)
          continue;
        if (!canPlace(present, newType))
          continue;
        StateType next = canonical(stateSet(cur, pos, newType));
        auto [it, inserted] = path.try_emplace(
            next, Transition{cur, (int8_t)pos, (int8_t)newType});
        if (inserted)
          q.push_back({next, curScore + (newType - curType)});
      }
    }
  }
  return path;
}

// Reverse BFS from a final state to base (all zeros).
// A tower can be removed if its placement requirements are met in the
// current state. The prior type at that position can be any other value,
// without needing to satisfy its own placement requirements.
static vector<pair<int, int>> buildPathFromReverseBfs(StateType finalState) {
  ankerl::unordered_dense::map<StateType, Transition> parentForward;
  ankerl::unordered_dense::map<StateType, Transition> parentReverse;
  parentForward.reserve(1ull << 26); // ~64M elements
  parentReverse.reserve(1ull << 26); // ~64M elements

  const size_t kPruneThreshold = 10000;
  const int kPruneLowDelta = 3;
  const int kPruneHighDelta = 3;

  auto countTowersLeft = [](StateType s) {
    int towersLeft = 0;
    for (int pos = 0; pos < g_numPos; pos++) {
      if (stateGet(s, pos) != 1)
        towersLeft++;
    }
    return towersLeft;
  };

  // Queue entries: (state, tower count) so pruning never calls countTowersLeft
  // over the full queue.
  using QueueEntry = std::pair<StateType, int>;

  // Prune: keep only frontier states that pass cutoff + their ancestors.
  // Build pruned map by walking backwards from kept frontier (no full parentMap
  // scan, no separate keepStates set).
  auto pruneSide = [&](auto &parentMap, std::deque<QueueEntry> &queue,
                       bool keepLowest) {
    if (parentMap.size() <= kPruneThreshold || queue.empty())
      return;

    int minTowers = g_numPos;
    int maxTowers = 0;
    for (const auto &entry : queue) {
      int tc = entry.second;
      minTowers = min(minTowers, tc);
      maxTowers = max(maxTowers, tc);
    }

    int cutoff = keepLowest ? (minTowers + kPruneLowDelta)
                            : (maxTowers - kPruneHighDelta);
    vector<QueueEntry> toProcess;
    toProcess.reserve(queue.size());
    for (const auto &entry : queue) {
      if ((keepLowest && entry.second <= cutoff) ||
          (!keepLowest && entry.second >= cutoff)) {
        toProcess.push_back(entry);
      }
    }

    if (toProcess.empty()) {
      auto bestIt = keepLowest ? min_element(queue.begin(), queue.end(),
                                             [](const auto &a, const auto &b) {
                                               return a.second < b.second;
                                             })
                               : max_element(queue.begin(), queue.end(),
                                             [](const auto &a, const auto &b) {
                                               return a.second > b.second;
                                             });
      toProcess.push_back(*bestIt);
    }

    using MapType = std::decay_t<decltype(parentMap)>;
    MapType pruned;
    pruned.reserve(toProcess.size() * 4);
    for (const auto &entry : toProcess)
      pruned[entry.first] = parentMap.at(entry.first);

    std::deque<QueueEntry> worklist(toProcess.begin(), toProcess.end());
    while (!worklist.empty()) {
      StateType cur = worklist.front().first;
      worklist.pop_front();
      auto it = parentMap.find(cur);
      if (it == parentMap.end() || it->second.getPos() == -1)
        continue;
      StateType parent = it->second.parent;
      if (parentMap.count(parent) && pruned.find(parent) == pruned.end()) {
        pruned[parent] = parentMap.find(parent)->second;
        worklist.push_back({parent, countTowersLeft(parent)});
      }
    }

    parentMap.swap(pruned);
    queue.assign(toProcess.begin(), toProcess.end());
  };

  std::deque<QueueEntry> qForward;
  std::deque<QueueEntry> qReverse;

  parentForward[g_baseState] = {g_baseState, -1, -1};
  qForward.push_back({g_baseState, 0});
  parentReverse[finalState] = {finalState, -1, -1};
  qReverse.push_back({finalState, countTowersLeft(finalState)});

  std::mutex coutMutex;
  std::atomic<bool> reverseFinished{false};
  size_t totalReverseVisited =
      0; // total states ever discovered by reverse (before pruning)

  // Result of expanding one state in reverse: (next state, transition, next
  // tower count).
  struct ReverseExpansion {
    StateType next;
    Transition trans;
    int nextTowers;
  };

  // Expand a range of (cur, curTowers) into a list of new states (no shared map
  // writes).
  auto reverseExpandChunk = [&](const QueueEntry *start, const QueueEntry *end,
                                vector<ReverseExpansion> *out) {
    for (; start != end; ++start) {
      StateType cur = start->first;
      int curTowers = start->second;
      if (cur == g_baseState)
        continue;
      for (int pos = 0; pos < g_numPos; pos++) {
        int curType = stateGet(cur, pos);
        if (curType == 1)
          continue;
        uint32_t present = neighborPresent(cur, pos);
        if (!canPlace(present, curType))
          continue;
        for (int prevType = 1; prevType <= g_maxType[pos]; prevType++) {
          if (prevType == curType)
            continue;
          StateType next = canonical(stateSet(cur, pos, prevType));
          int nextTowers =
              curTowers - (curType != 1 ? 1 : 0) + (prevType != 1 ? 1 : 0);
          out->push_back({next, Transition{cur, (int8_t)pos, (int8_t)curType},
                          nextTowers});
        }
      }
    }
  };

  // Reverse BFS: run to completion; each layer split between 2 threads.
  auto reverseTask = [&]() {
    size_t depth = 0;
    double lastReverseLayerMs = 0.0;
    while (!qReverse.empty()) {
      size_t layerCount = qReverse.size();
      {
        lock_guard<mutex> lock(coutMutex);
        cout << "[Reverse] Layer " << depth << ": expanding " << layerCount
             << " states ("
             << max(1, min(REVERSE_BFS_NUM_THREADS, (int)layerCount))
             << " threads)..." << endl;
      }
      auto layerStart = chrono::steady_clock::now();

      // Pop entire layer into a vector.
      vector<QueueEntry> layer;
      layer.reserve(layerCount);
      for (size_t i = 0; i < layerCount; i++) {
        layer.push_back(qReverse.front());
        qReverse.pop_front();
      }

      int minTowersLeft = g_numPos;
      int maxTowersLeft = 0;
      bool reachedBase = false;
      for (const auto &entry : layer) {
        if (entry.first == g_baseState) {
          reachedBase = true;
          break;
        }
        minTowersLeft = min(minTowersLeft, entry.second);
        maxTowersLeft = max(maxTowersLeft, entry.second);
      }
      if (reachedBase) {
        qReverse.clear();
        auto layerEnd = chrono::steady_clock::now();
        double expandMs =
            chrono::duration_cast<chrono::milliseconds>(layerEnd - layerStart)
                .count();
        lock_guard<mutex> lock(coutMutex);
        cout << "[Reverse] Layer " << depth << " done in " << expandMs
             << " ms (expand " << expandMs << " ms). Reached base state." << endl;
        break;
      }

      // Split layer across N threads; each expands its chunk into a local
      // vector.
      const int numThreads =
          max(1, min(REVERSE_BFS_NUM_THREADS, (int)layer.size()));
      const size_t chunkSize =
          (layer.size() + (size_t)numThreads - 1) / (size_t)numThreads;
      vector<vector<ReverseExpansion>> results((size_t)numThreads);
      for (int c = 0; c < numThreads; c++)
        results[c].reserve(layer.size() * 2);

      vector<thread> revThreads;
      for (int c = 0; c < numThreads; c++) {
        const QueueEntry *chunkStart = layer.data() + (size_t)c * chunkSize;
        const QueueEntry *chunkEnd =
            layer.data() + min((size_t)(c + 1) * chunkSize, layer.size());
        if (chunkStart < chunkEnd)
          revThreads.emplace_back(reverseExpandChunk, chunkStart, chunkEnd,
                                  &results[c]);
      }
      for (auto &t : revThreads)
        t.join();

      // Merge into shared map and queue (single-threaded, dedup via emplace).
      for (int c = 0; c < numThreads; c++) {
        for (const auto &r : results[c]) {
          if (parentReverse.emplace(r.next, r.trans).second) {
            totalReverseVisited++;
            qReverse.push_back({r.next, r.nextTowers});
          }
        }
      }

      auto layerEnd = chrono::steady_clock::now();
      double expandMs =
          chrono::duration_cast<chrono::milliseconds>(layerEnd - layerStart)
              .count();
      double pruneMs = 0.0;
      if (!qReverse.empty()) {
        auto pruneStart = chrono::steady_clock::now();
        pruneSide(parentReverse, qReverse, true);
        pruneMs = chrono::duration_cast<chrono::milliseconds>(
                      chrono::steady_clock::now() - pruneStart)
                      .count();
      }
      double totalMs = expandMs + pruneMs;
      lastReverseLayerMs = totalMs;
      {
        lock_guard<mutex> lock(coutMutex);
        cout << "[Reverse] Layer " << depth << " done in " << totalMs
             << " ms (expand " << expandMs << " ms, prune " << pruneMs
             << " ms). Queue size now " << qReverse.size() << ", map size "
             << parentReverse.size() << "." << endl;
      }
      depth++;
    }
    {
      lock_guard<mutex> lock(coutMutex);
      cout << "[Reverse] Queue empty. Reverse BFS finished." << endl;
    }
    reverseFinished = true;
  };

  // Forward BFS: stops as soon as reverse finishes (check inside layer so we
  // don't wait).
  auto forwardTask = [&]() {
    size_t forwardDepth = 0;
    while (!qForward.empty()) {
      if (reverseFinished) {
        lock_guard<mutex> lock(coutMutex);
        cout << "[Forward] Reverse finished; stopping immediately." << endl;
        break;
      }
      size_t forwardLayerCount = qForward.size();
      {
        lock_guard<mutex> lock(coutMutex);
        cout << "[Forward] Layer " << forwardDepth << ": expanding "
             << forwardLayerCount << " states..." << endl;
      }
      auto forwardLayerStart = chrono::steady_clock::now();
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
          minForwardTowersLeft = min(minForwardTowersLeft, curTowers);
          maxForwardTowersLeft = max(maxForwardTowersLeft, curTowers);
        }

        for (int pos = 0; pos < g_numPos; pos++) {
          int curType = stateGet(cur, pos);
          if (curType == 5)
            continue;
          uint32_t present = neighborPresent(cur, pos);
          for (int newType = 1; newType <= g_maxType[pos]; newType++) {
            if (newType == curType)
              continue;
            if (!canPlace(present, newType))
              continue;
            StateType next = canonical(stateSet(cur, pos, newType));
            if (parentForward
                    .emplace(next,
                             Transition{cur, (int8_t)pos, (int8_t)newType})
                    .second) {
              int nextTowers =
                  curTowers - (curType != 1 ? 1 : 0) + (newType != 1 ? 1 : 0);
              qForward.push_back({next, nextTowers});
            }
          }
        }
      }

      auto forwardLayerEnd = chrono::steady_clock::now();
      double expandMs = chrono::duration_cast<chrono::milliseconds>(
                            forwardLayerEnd - forwardLayerStart)
                            .count();
      double pruneMs = 0.0;
      if (!stoppedByReverse && !qForward.empty()) {
        auto pruneStart = chrono::steady_clock::now();
        pruneSide(parentForward, qForward, false);
        pruneMs = chrono::duration_cast<chrono::milliseconds>(
                      chrono::steady_clock::now() - pruneStart)
                      .count();
      }
      double totalMs = expandMs + pruneMs;
      {
        lock_guard<mutex> lock(coutMutex);
        cout << "[Forward] Layer " << forwardDepth << " done in " << totalMs
             << " ms (expand " << expandMs << " ms, prune " << pruneMs
             << " ms). Queue size now " << qForward.size() << ", map size "
             << parentForward.size();
        if (stoppedByReverse)
          cout << " (stopped mid-layer: reverse finished)";
        cout << "." << endl;
      }
      if (stoppedByReverse)
        break;
      forwardDepth++;
    }
    {
      lock_guard<mutex> lock(coutMutex);
      cout
          << "[Forward] Forward BFS finished (queue empty or reverse finished)."
          << endl;
    }
  };

  {
    lock_guard<mutex> lock(coutMutex);
    cout << "[Main] Starting reverse BFS (to completion) and forward BFS in "
            "parallel. Forward will stop as soon as reverse finishes."
         << endl;
  }
  thread tReverse(reverseTask);
  thread tForward(forwardTask);
  tReverse.join();
  {
    lock_guard<mutex> lock(coutMutex);
    cout << "[Main] Reverse finished. Not waiting for forward to complete; "
            "stopping forward and finding best intercept with current state."
         << endl;
  }
  tForward.join();
  {
    lock_guard<mutex> lock(coutMutex);
    cout << "[Main] Both threads done. Finding best intercept..." << endl;
  }

  // Best intercept: state that appears in both maps with minimum total path
  // length.
  auto forwardDepthOf = [&parentForward](StateType s) {
    int d = 0;
    StateType cur = s;
    while (parentForward.count(cur) && parentForward.at(cur).getPos() != -1) {
      cur = parentForward.at(cur).parent;
      d++;
    }
    return d;
  };
  auto reverseDepthOf = [&parentReverse, &finalState](StateType s) {
    int d = 0;
    StateType cur = s;
    while (parentReverse.count(cur) && parentReverse.at(cur).getPos() != -1) {
      cur = parentReverse.at(cur).parent;
      d++;
    }
    return d;
  };

  StateType meetState = g_baseState;
  int bestTotal = INT_MAX;
  size_t interceptCount = 0;
  for (const auto &kv : parentForward) {
    StateType s = kv.first;
    if (parentReverse.count(s) == 0)
      continue;
    interceptCount++;
    int fd = forwardDepthOf(s);
    int rd = reverseDepthOf(s);
    int total = fd + rd;
    if (total < bestTotal) {
      bestTotal = total;
      meetState = s;
    }
  }

  {
    lock_guard<mutex> lock(coutMutex);
    cout << "[Main] Found " << interceptCount
         << " states in both maps. Best intercept has forward depth "
         << forwardDepthOf(meetState) << ", reverse depth "
         << reverseDepthOf(meetState) << ", total " << bestTotal << " steps."
         << endl;
  }

  if (interceptCount == 0) {
    lock_guard<mutex> lock(coutMutex);
    cout << "[Main] No intercept: no state reached by both searches. Cannot "
            "reconstruct path."
         << endl;
    return {};
  }

  // Redirect output to file if baseline mode is active (after BFS table
  // completes)
  if (g_baselineFile && g_baselineFile->is_open()) {
    g_originalCoutBuf = cout.rdbuf();
    cout.rdbuf(g_baselineFile->rdbuf());
  }

  cout << "\nReverse BFS visited " << totalReverseVisited
       << " states (map size after pruning: " << parentReverse.size() << ").\n";

  vector<pair<int, int>> steps;
  vector<pair<int, int>> stepsForward;
  StateType cur = meetState;
  while (parentForward.at(cur).getPos() != -1) {
    const auto &tr = parentForward.at(cur);
    stepsForward.push_back({tr.getPos(), tr.getType()});
    cur = tr.parent;
  }
  reverse(stepsForward.begin(), stepsForward.end());

  cur = meetState;
  while (parentReverse.at(cur).getPos() != -1) {
    StateType parent = parentReverse.at(cur).parent;
    // canonical() may apply different symmetry transforms to cur and parent,
    // so their indices can differ at multiple positions. Search for the step
    // (pos, type) such that canonical(stateSet(cur, pos, type)) == parent.
    bool foundStep = false;
    for (int p = 0; p < g_numPos && !foundStep; p++) {
      for (int t = 1; t <= g_maxType[p] && !foundStep; t++) {
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

  steps.insert(steps.begin(), stepsForward.begin(), stepsForward.end());
  return steps;
}

// ── Sliding Window Optimization ─────────────────────────────────────────────-
// Processes grid using 3x3 windows with stride=1 along a JPEG-like zigzag path
// Each window respects constraints from the entire grid (in+out of window)

// Build zigzag path starting from a given corner
// corner: 0=top-left, 1=top-right, 2=bottom-left, 3=bottom-right
static vector<int> buildZigzagPath(int corner = 0) {
  int R = GRID.rows, C = GRID.cols;
  vector<vector<pair<int, int>>> diagonals;

  // Group by diagonal - direction depends on corner
  for (int r = 0; r < R; r++) {
    for (int c = 0; c < C; c++) {
      pair<int, int> coord = {r, c};
      if (!GRID.coordToIdx.count(coord))
        continue;

      int diag;
      switch (corner) {
      case 0:
        diag = r + c;
        break; // top-left
      case 1:
        diag = r + (C - 1 - c);
        break; // top-right
      case 2:
        diag = (R - 1 - r) + c;
        break; // bottom-left
      case 3:
        diag = (R - 1 - r) + (C - 1 - c);
        break; // bottom-right
      default:
        diag = r + c;
        break;
      }
      if (diag >= (int)diagonals.size())
        diagonals.resize(diag + 1);
      diagonals[diag].push_back({r, c});
    }
  }

  vector<int> path;
  for (int d = 0; d < (int)diagonals.size(); d++) {
    auto &diag = diagonals[d];
    // Alternate direction for each diagonal (zigzag effect)
    if (d % 2 == 1)
      reverse(diag.begin(), diag.end());
    for (auto [r, c] : diag) {
      path.push_back(GRID.coordToIdx[{r, c}]);
    }
  }
  return path;
}

// BFS within a window: exhaustively explore ALL reachable states
// by modifying only positions within the window
static void exploreWindow(StateType startState, int windowCenterIdx,
                          vector<StateType> &newStates, int &bestScore,
                          StateType &bestState, long long &explored) {
  // Get kernel pattern window around position
  auto [cr, cc] = GRID.positions[windowCenterIdx];

  // Parse kernel pattern to get offsets
  int patternRows = (int)WINDOW_KERNEL_PATTERN.size();
  int patternCols = 0;
  for (const auto &row : WINDOW_KERNEL_PATTERN) {
    patternCols = max(patternCols, (int)row.size());
  }
  int centerR = patternRows / 2;
  int centerC = patternCols / 2;

  vector<int> windowPositions;
  for (int pr = 0; pr < patternRows; pr++) {
    for (int pc = 0; pc < (int)WINDOW_KERNEL_PATTERN[pr].size(); pc++) {
      if (WINDOW_KERNEL_PATTERN[pr][pc] == '.') {
        int dr = pr - centerR;
        int dc = pc - centerC;
        auto it = GRID.coordToIdx.find({cr + dr, cc + dc});
        if (it != GRID.coordToIdx.end()) {
          windowPositions.push_back(it->second);
        }
      }
    }
  }

  // Full BFS within window - explore ALL states reachable by modifying window
  // positions
  robin_hood::unordered_set<StateType> localVisited;
  localVisited.reserve(1ull << 24); // ~16M elements per window
  struct Entry {
    StateType state;
    int score;
    int depth;
  };
  vector<Entry> q;

  localVisited.insert(startState);
  q.push_back({startState, scoreGrid(startState), 0});
  size_t head = 0;

  while (head < q.size()) {
    auto [cur, curScore, depth] = q[head++];
    explored++;

    // Stop expanding this branch if we've reached max depth
    if (depth >= MAX_DEPTH)
      continue;

    if (curScore > bestScore) {
      bestScore = curScore;
      bestState = cur;
      cout << "New best: " << bestScore << "\n";
      printGrid(bestState);
    }

    // Try modifying each position in the window
    for (int pos : windowPositions) {
      if (pos >= g_numPos)
        continue;
      int curType = stateGet(cur, pos);

      // Try all valid types for this position
      uint32_t present = neighborPresent(cur, pos);
      for (int newType = 1; newType <= g_maxType[pos]; newType++) {
        if (newType == curType)
          continue;
        if (!canPlace(present, newType))
          continue;

        StateType next = stateSet(cur, pos, newType);

        // Only track within this window to prevent infinite loops
        if (localVisited.insert(next).second) {
          int nextScore = curScore + (newType - curType);
          q.push_back({next, nextScore, depth + 1});
          newStates.push_back(next);
        }
      }
    }
  }
  // localVisited is destroyed here - no global state persists
}

static SearchResult runSlidingWindow() {
  SearchResult result{g_baseState, 0, 0, 0};
  result.bestState = g_baseState;
  result.bestScore = 0;
  result.explored = 0;

  cout << "Sliding window - sequential greedy approach:\n";
  for (const auto &row : WINDOW_KERNEL_PATTERN) {
    cout << "  " << row << "\n";
  }

  const char *cornerNames[] = {"top-left", "top-right", "bottom-left",
                               "bottom-right"};

  struct WindowStats {
    int posIdx;
    int row, col;
    int kernelSize;
    long long statesExplored;
    long long timeMs;
    bool improved;
    int scoreAfter;
  };

  StateType currentState = g_baseState;
  int currentScore = 0;
  int roundNum = 0;

  while (true) {
    roundNum++;
    int scoreAtRoundStart = currentScore;

    cout << "============ Round " << roundNum
         << " (starting score: " << currentScore << ") ============\n\n";

    int corner = 1; // top-right
    vector<int> zigzag = buildZigzagPath(corner);

    int bestScore = currentScore;
    StateType bestState = currentState;
    long long explored = 0;

    cout << "=== Pass " << (char)('A' + corner) << ": " << cornerNames[corner]
         << " (path length: " << zigzag.size() << ") ===\n";

    auto passStart = chrono::steady_clock::now();

    vector<WindowStats> stats;

    for (int windowIdx = 0; windowIdx < (int)zigzag.size(); windowIdx++) {
      vector<StateType> exploredStates;
      long long exploredBefore = explored;

      int posIdx = zigzag[windowIdx];
      auto [cr, cc] = GRID.positions[posIdx];

      auto t0 = chrono::steady_clock::now();

      exploreWindow(currentState, posIdx, exploredStates, bestScore, bestState,
                    explored);

      auto t1 = chrono::steady_clock::now();
      auto ms = chrono::duration_cast<chrono::milliseconds>(t1 - t0).count();

      bool improved = (bestState != currentState);
      if (improved) {
        currentState = bestState;
      }

      stats.push_back({posIdx, cr, cc, (int)exploredStates.size(),
                       explored - exploredBefore, ms, improved, bestScore});
    }

    currentScore = bestScore;

    auto passEnd = chrono::steady_clock::now();
    auto passMs =
        chrono::duration_cast<chrono::milliseconds>(passEnd - passStart)
            .count();

    // Print summary table
    cout << "\n┌─────┬─────────┬──────────────┬────────────────┬──────────┬────"
            "──┬───────"
            "┐\n";
    cout << "│ Pos │  (r,c)  │    Kernel    │     States     │   Time   │ Imp  "
            "│ "
            "Score │\n";
    cout << "├─────┼─────────┼──────────────┼────────────────┼──────────┼──────"
            "┼──────"
            "─┤\n";
    for (const auto &s : stats) {
      printf("│ %3d │ (%d, %d)  │ %12d │ %14lld │ %6lldms │  %s  │ %5d │\n",
             s.posIdx, s.row, s.col, s.kernelSize, s.statesExplored, s.timeMs,
             s.improved ? "Y" : "N", s.scoreAfter);
    }
    cout << "└─────┴─────────┴──────────────┴────────────────┴──────────┴──────"
            "┴──────"
            "─┘\n";
    cout << "Pass " << (char)('A' + corner) << " (" << cornerNames[corner]
         << "): score " << currentScore << " (" << passMs << "ms, " << explored
         << " states)\n\n";

    result.explored += explored;

    cout << "Round " << roundNum << " complete: score " << scoreAtRoundStart
         << " -> " << currentScore << "\n";

    if (currentScore <= scoreAtRoundStart) {
      cout << "No improvement in round " << roundNum << ", stopping.\n\n";
      break;
    } else {
      cout << "Improved by " << (currentScore - scoreAtRoundStart)
           << ", starting another round.\n\n";
    }
  }

  result.bestScore = currentScore;
  result.bestState = currentState;

  cout << "\nFinal score after " << roundNum
       << " round(s): " << result.bestScore << "\n";

  // ── Build path reconstruction: reverse BFS ─────────────────────────────
  StateType finalState = canonical(result.bestState);
  vector<pair<int, int>> buildPath = buildPathFromReverseBfs(finalState);
  if (buildPath.empty()) {
    cout << "\nBuild path not found by reverse BFS.\n";
    cout << "\nFinal score: " << result.bestScore << "\n";
    result.uniqueStates = 0;
    return result;
  }

  printExpandedPath(expandPathForEmptyBoard(buildPath));

  cout << "\nFinal score: " << result.bestScore << "\n";
  result.uniqueStates = 0;

  // Restore cout if it was redirected
  if (g_originalCoutBuf) {
    cout.rdbuf(g_originalCoutBuf);
    g_originalCoutBuf = nullptr;
  }

  return result;
}

static SearchResult runIDDFS() {
  SearchResult result{0, 0, 0, 1};
  robin_hood::unordered_set<StateType> globalSeen;
  globalSeen.reserve(1ull << 24); // ~16M elements
  globalSeen.insert(g_baseState);

  for (int depthLimit = 1;; depthLimit++) {
    robin_hood::unordered_set<StateType> seen;
    seen.reserve(1ull << 20); // ~1M elements per depth
    struct Entry {
      StateType state;
      int score;
      int depth;
    };
    vector<Entry> stack;
    size_t newStates = 0;

    stack.push_back({g_baseState, 0, 0});
    seen.insert(g_baseState);

    while (!stack.empty()) {
      auto cur = stack.back();
      stack.pop_back();
      result.explored++;

      if (cur.score > result.bestScore) {
        result.bestScore = cur.score;
        result.bestState = cur.state;
        cout << "New best score: " << result.bestScore << "\n";
        printGrid(cur.state);
      }

      if (cur.depth >= depthLimit)
        continue;

      for (int pos = 0; pos < g_numPos; pos++) {
        int curType = stateGet(cur.state, pos);
        if (curType == 5)
          continue;
        uint32_t present = neighborPresent(cur.state, pos);
        for (int newType = 1; newType <= g_maxType[pos]; newType++) {
          if (newType == curType)
            continue;
          if (!canPlace(present, newType))
            continue;

          StateType next = canonical(stateSet(cur.state, pos, newType));
          int nextScore = cur.score + (newType - curType);
          if (seen.insert(next).second) {
            if (globalSeen.insert(next).second)
              newStates++;
            stack.push_back({next, nextScore, cur.depth + 1});
          }
        }
      }
    }

    result.uniqueStates = globalSeen.size();
    cout << "Depth " << depthLimit << " complete. New states: " << newStates
         << "\n";
    if (newStates == 0)
      break;
  }

  return result;
}

// ── Public entrypoint for a single level/unlocked-type configuration ────────

void runSingleLevel() {
  cout << "Tower Blox Deluxe Solver\n";
  cout << "========================\n";
  cout << "Level: " << getCurrentLevel() << "\n";
  cout << "Tower types unlocked (1-" << MAX_UNLOCKED_TYPE << "):\n";
  for (int i = 1; i <= MAX_UNLOCKED_TYPE; i++)
    cout << "  " << COLORS[i] << BLACK_TEXT << i << RESET << " ("
         << COLOR_NAMES[i] << ")\n";
  if (SEARCH_MODE == 0)
    cout << "\nRunning BFS over all reachable states...\n\n";
  else if (SEARCH_MODE == 1)
    cout << "\nRunning DFS over all reachable states...\n\n";
  else if (SEARCH_MODE == 2)
    cout << "\nRunning IDDFS over all reachable states...\n\n";
  else
    cout << "\nRunning sliding window optimization (sequential greedy)...\n\n";

  buildLookupTables();
  buildSymmetries();

  SearchResult result{0, 0, 0, 0};
  robin_hood::unordered_map<StateType, Transition> visited;
  visited.reserve(1ull << 24); // ~16M elements

  if (SEARCH_MODE == 3) {
    result = runSlidingWindow();
  } else if (SEARCH_MODE == 2) {
    result = runIDDFS();
  } else {
    // BFS/DFS: every state is visited exactly once.
    // For BFS, the vector is used as a queue with a head index.
    // For DFS, the vector is used as a stack via push_back/pop_back.
    struct Entry {
      StateType state;
      int score;
    };
    vector<Entry> q;
    size_t head = 0;

    StateType startState = g_baseState;
    visited[startState] = {g_baseState, -1, -1};
    q.push_back({startState, 0});

    result.bestState = startState;
    result.bestScore = 0;
    result.explored = 0;

    while ((SEARCH_MODE == 0) ? (head < q.size()) : !q.empty()) {
      Entry curEntry;
      if (SEARCH_MODE == 0) {
        curEntry = q[head++];
      } else {
        curEntry = q.back();
        q.pop_back();
      }
      auto [cur, curScore] = curEntry;
      result.explored++;

      if (curScore > result.bestScore) {
        result.bestScore = curScore;
        result.bestState = cur;
        cout << "New best score: " << result.bestScore << "\n";
        printGrid(cur);
      }

      // Try every possible single-position change
      for (int pos = 0; pos < g_numPos; pos++) {
        int curType = stateGet(cur, pos);
        if (curType == 5)
          continue;
        uint32_t present = neighborPresent(cur, pos); // computed once per pos
        for (int newType = 1; newType <= g_maxType[pos]; newType++) {
          if (newType == curType)
            continue;
          if (!canPlace(present, newType))
            continue;

          StateType next = canonical(stateSet(cur, pos, newType));
          int nextScore = curScore + (newType - curType);
          auto [it, inserted] = visited.try_emplace(
              next, Transition{cur, (int8_t)pos, (int8_t)newType});
          if (inserted) {
            q.push_back({next, nextScore});
          }
        }
      }

      if (result.explored % 2000000 == 0)
        cout << "Explored: " << result.explored
             << " | Visited: " << visited.size()
             << " | Best: " << result.bestScore << "\r" << flush;
    }

    result.uniqueStates = visited.size();
  }

  if (SEARCH_MODE == 3)
    cout << "\nSliding window complete.\n";
  else if (SEARCH_MODE == 2)
    cout << "\nIDDFS complete.\n";
  else if (SEARCH_MODE == 1)
    cout << "\nDFS complete.\n";
  else
    cout << "\nBFS complete.\n";
  cout << "States explored : " << result.explored << "\n";
  cout << "Unique states   : " << result.uniqueStates << "\n";
  cout << "Best score      : " << result.bestScore << "\n\n";
  cout << "Optimal layout:\n";
  printGrid(result.bestState);

  // Path reconstruction only available for BFS
  if (SEARCH_MODE == 0) {
    printBuildPath(result.bestState, visited);
  } else if (SEARCH_MODE == 3) {
    // Path already printed by runSlidingWindow
  } else {
    cout << "\nFinding shortest build path via BFS...\n";
    auto pathMap = bfsShortestPath(result.bestState);
    printBuildPath(result.bestState, pathMap);
  }
}
