cmake_minimum_required(VERSION 3.10)
project(TowerBloxSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O3 /GL /arch:AVX2")
  set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")
  add_compile_options(/MP)
else()
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")
endif()

# Ankerl unordered_dense library
include(FetchContent)
FetchContent_Declare(unordered_dense
  GIT_REPOSITORY https://github.com/martinus/unordered_dense.git
  GIT_TAG v4.4.0
)
FetchContent_MakeAvailable(unordered_dense)

add_executable(tower_solver main.cpp)
target_include_directories(tower_solver PRIVATE ${unordered_dense_SOURCE_DIR}/include)
