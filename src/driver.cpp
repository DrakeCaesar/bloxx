#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

#ifdef _WIN32
#include <windows.h>
#endif

#include "config.h"
#include "display.h"
#include "grid_config.h"
#include "search.h"

using namespace std;

int main() {
#ifdef _WIN32
  // Set console to UTF-8 mode for proper Unicode display on Windows
  SetConsoleOutputCP(CP_UTF8);
#endif

  // Process configured level range and all unlocked tower types

  const int numLevels = LEVEL_END - LEVEL_START + 1;
  const int numTypes = UNLOCKED_TYPE_END - UNLOCKED_TYPE_START + 1;
  const int totalConfigs = numLevels * numTypes;
  int currentIndex = 0;
  auto totalStart = chrono::steady_clock::now();

  for (int level = LEVEL_START; level <= LEVEL_END; level++) {
    for (int unlockedType = UNLOCKED_TYPE_START; unlockedType <= UNLOCKED_TYPE_END; unlockedType++) {
      currentIndex++;
      cout << "Config " << currentIndex << "/" << totalConfigs << ": level " << level << " unlocked 1-" << unlockedType
           << "  " << progressBar((size_t)currentIndex, (size_t)totalConfigs) << "\n";

      // If running from build/, write baselines one level up; otherwise current
      // dir
      string prefix = "";
      if (std::filesystem::current_path().filename() == "build")
        prefix = "../";
      string filename = prefix + "level" + to_string(level) + "_type" + to_string(unlockedType) + "_baseline.txt";
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

      auto runStart = chrono::steady_clock::now();
      // Run the solver (console output for sliding window, file output for
      // reverse BFS)
      runSingleLevel();
      auto runEnd = chrono::steady_clock::now();
      auto runMs = chrono::duration_cast<chrono::milliseconds>(runEnd - runStart).count();

      // Clear global file pointer
      g_baselineFile = nullptr;
      outFile.close();

      cout << "  Level " << level << " (unlocked: 1-" << unlockedType << ") complete in "
           << formatTimeMs((long long)runMs) << " -> " << filename << "\n";
    }
  }

  auto totalEnd = chrono::steady_clock::now();
  auto totalMs = chrono::duration_cast<chrono::milliseconds>(totalEnd - totalStart).count();
  cout << "All " << totalConfigs << " baseline files generated in " << formatTimeMs((long long)totalMs) << ".\n";

  return 0;
}
