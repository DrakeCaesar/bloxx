#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
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
// Layout:
//   X 0 1 X
//   2 3 4 5
//   6 7 8 9
//   X A B X
// Total 12 positions

struct Grid {
  int towers[12];

  Grid() { fill(towers, towers + 12, 0); }

  int &operator[](int i) { return towers[i]; }
  const int &operator[](int i) const { return towers[i]; }

  int getScore() const {
    int score = 0;
    for (int i = 0; i < 12; i++) {
      score += towers[i];
    }
    return score;
  }
};

// Get neighbors for each position
vector<int> getNeighbors(int pos) {
  // Position mapping:
  //   X  0  1  X
  //   2  3  4  5
  //   6  7  8  9
  //   X  10 11 X

  vector<int> neighbors;

  switch (pos) {
  case 0:
    neighbors = {1, 2, 3};
    break;
  case 1:
    neighbors = {0, 3, 4, 5};
    break;
  case 2:
    neighbors = {0, 3, 6};
    break;
  case 3:
    neighbors = {0, 1, 2, 4, 6, 7};
    break;
  case 4:
    neighbors = {1, 3, 5, 7, 8};
    break;
  case 5:
    neighbors = {1, 4, 8, 9};
    break;
  case 6:
    neighbors = {2, 3, 7, 10};
    break;
  case 7:
    neighbors = {3, 4, 6, 8, 10, 11};
    break;
  case 8:
    neighbors = {4, 5, 7, 9, 11};
    break;
  case 9:
    neighbors = {5, 8, 11};
    break;
  case 10:
    neighbors = {6, 7, 11};
    break;
  case 11:
    neighbors = {7, 8, 9, 10};
    break;
  }

  return neighbors;
}

// Check if a tower placement is valid
bool isValid(const Grid &grid, int pos) {
  int towerType = grid[pos];

  // Type 0 has no requirements
  if (towerType == 0)
    return true;

  vector<int> neighbors = getNeighbors(pos);
  set<int> requiredTypes;

  // Build set of required neighboring tower types
  for (int i = 0; i < towerType; i++) {
    requiredTypes.insert(i);
  }

  // Check if all required types are present in neighbors
  set<int> foundTypes;
  for (int n : neighbors) {
    foundTypes.insert(grid[n]);
  }

  for (int req : requiredTypes) {
    if (foundTypes.find(req) == foundTypes.end()) {
      return false;
    }
  }

  return true;
}

// Check if entire grid is valid
bool isValidGrid(const Grid &grid) {
  for (int i = 0; i < 12; i++) {
    if (!isValid(grid, i)) {
      return false;
    }
  }
  return true;
}

// Print the grid
void printGrid(const Grid &grid, int score) {
  cout << "\n=== Layout Found! Score: " << score << " ===\n";
  cout << "  X  " << COLORS[grid[0]] << grid[0] << RESET << "  "
       << COLORS[grid[1]] << grid[1] << RESET << "  X\n";
  cout << "  " << COLORS[grid[2]] << grid[2] << RESET << "  " << COLORS[grid[3]]
       << grid[3] << RESET << "  " << COLORS[grid[4]] << grid[4] << RESET
       << "  " << COLORS[grid[5]] << grid[5] << RESET << "\n";
  cout << "  " << COLORS[grid[6]] << grid[6] << RESET << "  " << COLORS[grid[7]]
       << grid[7] << RESET << "  " << COLORS[grid[8]] << grid[8] << RESET
       << "  " << COLORS[grid[9]] << grid[9] << RESET << "\n";
  cout << "  X  " << COLORS[grid[10]] << grid[10] << RESET << "  "
       << COLORS[grid[11]] << grid[11] << RESET << "  X\n";
  cout << "================================\n";
}

// Brute force all possibilities
void bruteForce() {
  Grid grid;
  int maxScore = -1;
  long long totalConfigs = 0;
  long long validConfigs = 0;

  // Try all combinations: 5^12 possibilities
  long long maxCombinations = 1;
  for (int i = 0; i < 12; i++)
    maxCombinations *= 5;

  cout << "Starting brute force search...\n";
  cout << "Total combinations to check: " << maxCombinations << "\n\n";

  for (long long combo = 0; combo < maxCombinations; combo++) {
    totalConfigs++;

    // Convert combo number to grid configuration
    long long temp = combo;
    for (int i = 0; i < 12; i++) {
      grid[i] = temp % 5;
      temp /= 5;
    }

    // Check if valid
    if (isValidGrid(grid)) {
      validConfigs++;
      int score = grid.getScore();

      // Print if better score found
      if (score > maxScore) {
        maxScore = score;
        printGrid(grid, score);
      }
    }

    // Progress indicator
    if (combo % 1000000 == 0) {
      cout << "Progress: " << (combo * 100 / maxCombinations) << "%\r" << flush;
    }
  }

  cout << "\n\nSearch complete!\n";
  cout << "Total configurations checked: " << totalConfigs << "\n";
  cout << "Valid configurations found: " << validConfigs << "\n";
  cout << "Best score: " << maxScore << "\n";
}

int main() {
  cout << "Tower Blox Deluxe City Layout Solver\n";
  cout << "=====================================\n";
  cout << "Tower Types:\n";
  cout << "  " << COLORS[0] << "0" << RESET << " (Blue)   - No requirements\n";
  cout << "  " << COLORS[1] << "1" << RESET << " (Red)    - Needs 0 neighbor\n";
  cout << "  " << COLORS[2] << "2" << RESET
       << " (Green)  - Needs 0 and 1 neighbors\n";
  cout << "  " << COLORS[3] << "3" << RESET
       << " (Yellow) - Needs 0, 1, and 2 neighbors\n";
  cout << "  " << COLORS[4] << "4" << RESET
       << " (White)  - Needs 0, 1, 2, and 3 neighbors\n";
  cout << "\n";

  bruteForce();

  return 0;
}
