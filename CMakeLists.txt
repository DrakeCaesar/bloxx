cmake_minimum_required(VERSION 3.10)
project(TowerBloxSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O3 /GL /arch:AVX2")
  set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")
  add_compile_options(/MP)
  # RelWithDebInfo: use for profiling in VS (full PDB, walkable stacks)
  # Build: cmake -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo && cmake --build build --config RelWithDebInfo
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} /O2 /Oi /arch:AVX2 /Oy-")
  set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO} /DEBUG:FULL")
else()
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -flto")
  set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -flto")
endif()

# Robin-hood hashing library
include(FetchContent)
FetchContent_Declare(robin_hood
  GIT_REPOSITORY https://github.com/martinus/robin-hood-hashing.git
  GIT_TAG 3.11.5
)
FetchContent_MakeAvailable(robin_hood)

# Ankerl unordered_dense library
FetchContent_Declare(unordered_dense
  GIT_REPOSITORY https://github.com/martinus/unordered_dense.git
  GIT_TAG v4.8.1
)
FetchContent_MakeAvailable(unordered_dense)

add_executable(tower_solver src/driver.cpp src/search.cpp)
target_link_libraries(tower_solver PRIVATE robin_hood::robin_hood)
target_include_directories(tower_solver PRIVATE ${unordered_dense_SOURCE_DIR}/include)
