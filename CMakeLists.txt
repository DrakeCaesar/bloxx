cmake_minimum_required(VERSION 3.10)
project(TowerBloxSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /Oi /Ot /GL /fp:fast /GS- /Gy")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG /OPT:REF /OPT:ICF")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /arch:AVX2")

# Parallel compilation
add_compile_options(/MP)

# Robin-hood hashing library
include(FetchContent)
FetchContent_Declare(robin_hood
  GIT_REPOSITORY https://github.com/martinus/robin-hood-hashing.git
  GIT_TAG 3.11.5
)
FetchContent_MakeAvailable(robin_hood)

add_executable(tower_solver main.cpp)
target_link_libraries(tower_solver PRIVATE robin_hood::robin_hood)
