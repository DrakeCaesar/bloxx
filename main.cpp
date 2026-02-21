#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <vector>

using namespace std;

// ANSI color codes
const string COLORS[] = {
    "\033[34m", // Blue for type 0
    "\033[31m", // Red for type 1
    "\033[32m", // Green for type 2
    "\033[33m", // Yellow for type 3
    "\033[37m"  // White for type 4
};
const string RESET = "\033[0m";

// Grid is 4x4 without corners
//   X  0  1  X
//   2  3  4  5
//   6  7  8  9
//   X  10 11 X
// Total 12 positions (all start as type 0)

struct State {
  int towers[12];

  State() {
    for (int i = 0; i < 12; i++)
      towers[i] = 0; // All start as type 0 (base state)
  }

  int &operator[](int i) { return towers[i]; }
  const int &operator[](int i) const { return towers[i]; }

  int getScore() const {
    int score = 0;
    for (int i = 0; i < 12; i++) {
      score += towers[i];
    }
    return score;
  }

  bool operator<(const State &other) const {
    for (int i = 0; i < 12; i++) {
      if (towers[i] != other[i])
        return towers[i] < other[i];
    }
    return false;
  }
};

// Neighbor lookup
const vector<int> NEIGHBORS[12] = {
    {1, 2, 3},            // 0
    {0, 3, 4, 5},         // 1
    {0, 3, 6},            // 2
    {0, 1, 2, 4, 6, 7},   // 3
    {1, 3, 5, 7, 8},      // 4
    {1, 4, 8, 9},         // 5
    {2, 3, 7, 10},        // 6
    {3, 4, 6, 8, 10, 11}, // 7
    {4, 5, 7, 9, 11},     // 8
    {5, 8, 11},           // 9
    {6, 7, 11},           // 10
    {7, 8, 9, 10}         // 11
};

// Check if we can place a tower at position
bool canPlace(const State &state, int pos, int type) {
  if (type == 0)
    return true;

  // Need neighbors of types 0, 1, ..., type-1
  bool hasType[5] = {false};
  for (int neighbor : NEIGHBORS[pos]) {
    hasType[state[neighbor]] = true;
  }

  // Check all required types are present
  for (int i = 0; i < type; i++) {
    if (!hasType[i])
      return false;
  }

  return true;
}

// Print state
void printState(const State &state, int score) {
  cout << "\n=== New Best Score: " << score << " ===\n";
  cout << "  X  ";
  cout << COLORS[state[0]] << state[0] << RESET << "  ";
  cout << COLORS[state[1]] << state[1] << RESET << "  X\n";

  cout << "  ";
  for (int i = 2; i <= 5; i++) {
    cout << COLORS[state[i]] << state[i] << RESET << "  ";
  }
  cout << "\n";

  cout << "  ";
  for (int i = 6; i <= 9; i++) {
    cout << COLORS[state[i]] << state[i] << RESET << "  ";
  }
  cout << "\n";

  cout << "  X  ";
  cout << COLORS[state[10]] << state[10] << RESET << "  ";
  cout << COLORS[state[11]] << state[11] << RESET << "  X\n";
  cout << "=============================\n";
}

// Global tracking
int bestScore = -1;
State bestState;
long long statesExplored = 0;
const int MAX_DEPTH = 10; // Depth limit instead of visited set

// Simple DFS exploration with depth limit
void dfs(State &current, int depth) {
  statesExplored++;

  // Check if this is a new best
  int score = current.getScore();
  if (score > bestScore) {
    bestScore = score;
    bestState = current;
    printState(current, score);
  }

  // Progress indicator
  if (statesExplored % 100000 == 0) {
    cout << "Explored: " << statesExplored << " states | Best: " << bestScore
         << " | Depth: " << depth << "\r" << flush;
  }

  // Stop if depth limit reached
  if (depth >= MAX_DEPTH)
    return;

  // Try upgrading or downgrading each position
  for (int pos = 0; pos < 12; pos++) {
    int currentType = current[pos];

    // Try all different types (0-4)
    for (int newType = 0; newType <= 4; newType++) {
      if (newType != currentType && canPlace(current, pos, newType)) {
        current[pos] = newType;
        dfs(current, depth + 1);
        current[pos] = currentType; // Backtrack
      }
    }
  }
}

int main() {
  cout << "Tower Blox Deluxe Solver - Simple DFS\n";
  cout << "======================================\n";
  cout << "Tower Types:\n";
  cout << "  " << COLORS[0] << "0" << RESET
       << " (Blue)   - Base state (no requirements)\n";
  cout << "  " << COLORS[1] << "1" << RESET << " (Red)    - Needs 0 neighbor\n";
  cout << "  " << COLORS[2] << "2" << RESET
       << " (Green)  - Needs 0 and 1 neighbors\n";
  cout << "  " << COLORS[3] << "3" << RESET
       << " (Yellow) - Needs 0, 1, and 2 neighbors\n";
  cout << "  " << COLORS[4] << "4" << RESET
       << " (White)  - Needs 0, 1, 2, and 3 neighbors\n";
  cout << "\nNote: All positions start at type 0 (base state).\n";
  cout << "Algorithm tries upgrading/downgrading each position.\n";
  cout << "\nStarting DFS exploration from all type 0 (base state)...\n\n";

  State start;
  dfs(start, 0);

  cout << "\n\nExploration Complete!\n";
  cout << "States explored: " << statesExplored << "\n";
  cout << "Best score: " << bestScore << "\n\n";

  cout << "Final best layout:\n";
  printState(bestState, bestScore);

  return 0;
}
