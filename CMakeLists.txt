
cmake_minimum_required(VERSION 3.0)

project(MicroPP VERSION 0.1
  LANGUAGES C CXX Fortran)

# Use C++11
if (CMAKE_VERSION VERSION_LESS "3.1")
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
  endif ()
else ()
  set (CMAKE_CXX_STANDARD 11)
endif ()

option(TIMER "Enable time benchmarks instrumentation" OFF)

if (TIMER)
  add_definitions(-DTIMER)
endif()

# Include Directories (for all targets)
include_directories(include ${CMAKE_BINARY_DIR})

# wildcard all the sources in src
file(GLOB SOURCESLIB src/*.cpp src/*.f95)

add_library(micropp ${SOURCESLIB})

enable_testing()
add_subdirectory(test)
