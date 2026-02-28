#pragma once

#include <fstream>
#include <streambuf>

// High-level solver entrypoint for a single level and unlocked-type
void runSingleLevel();

// Baseline output redirection hooks used by the solver
extern std::ofstream *g_baselineFile;
extern std::streambuf *g_originalCoutBuf;

