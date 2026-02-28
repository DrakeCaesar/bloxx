#include <fstream>
#include <iostream>

#ifdef _WIN32
#include <windows.h>
#endif

#include "config.h"
#include "grid_config.h"
#include "search.h"

using namespace std;

int main() {
#ifdef _WIN32
  // Set console to UTF-8 mode for proper Unicode display on Windows
  SetConsoleOutputCP(CP_UTF8);
  setvbuf(stdout, nullptr, _IOFBF, 1000);
#endif

  // Process configured level range and all unlocked tower types
  cout << "═══════════════════════════════════════════════════════════════\n";
  cout << "       Tower Blox Deluxe Solver - Multi-Level Baseline\n";
  cout << "═══════════════════════════════════════════════════════════════\n\n";

  for (int level = LEVEL_START; level <= LEVEL_END; level++) {
    for (int unlockedType = UNLOCKED_TYPE_START;
         unlockedType <= UNLOCKED_TYPE_END; unlockedType++) {
      // Generate filename with both level and unlocked type
      string filename = "../level" + to_string(level) + "_type" +
                        to_string(unlockedType) + "_baseline.txt";
      ofstream outFile(filename, std::ios::binary);
      if (!outFile) {
        cerr << "Failed to open " << filename << "\n";
        continue;
      }

      // Set global file pointer (will be used during reverse BFS)
      g_baselineFile = &outFile;

      // Reinitialize for this level and tower type
      reinitializeGrid(level);
      GRID.computeMaxTypes(unlockedType);

      // Run the solver (console output for sliding window, file output for
      // reverse BFS)
      runSingleLevel();

      // Clear global file pointer
      g_baselineFile = nullptr;
      outFile.close();

      // Print progress to console
      cout << "✓ Level " << level << " (unlocked: 0-" << unlockedType
           << ") complete → " << filename << "\n";
    }
  }

  cout << "\n═══════════════════════════════════════════════════════════════\n";
  cout << "All 25 baseline files generated successfully!\n";
  cout << "═══════════════════════════════════════════════════════════════\n";

  return 0;
}

