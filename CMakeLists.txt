

cmake_minimum_required(VERSION 3.0)

project(MicroCpp VERSION 0.1
  DESCRIPTION "Description here!!"
  LANGUAGES C CXX Fortran)

# Use C++11
set(CMAKE_CXX_STANDARD 11)

# Include Directories (for all targets)
include_directories(include)

# wildcard all the sources in src
file(GLOB SOURCESLIB src/*.cpp)

# Library
add_library(microcpp ${SOURCESLIB})

# Enable auto create tests
enable_testing()
add_subdirectory(test)
