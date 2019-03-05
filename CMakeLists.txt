
cmake_minimum_required(VERSION 3.0)

project(MicroPP VERSION 0.1 LANGUAGES C CXX Fortran)

# Use C++11
if (CMAKE_VERSION VERSION_LESS "3.1")
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
  endif ()
else ()
  set (CMAKE_CXX_STANDARD 11)
endif ()

option(TIMER "Enable time benchmarks instrumentation" OFF)
option(CGDEBUG "Enable CG algorithm debug" OFF)
option(FILTER "Enable Filtering on Average Stress and Ctan" OFF)
option(OPENMP "Enables OpenMP paralelization" OFF)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif(NOT CMAKE_BUILD_TYPE)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")

if(OPENMP)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

if (TIMER)
	add_definitions(-DTIMER)
endif()

if (CGDEBUG)
	add_definitions(-DCGDEBUG)
endif()

if (FILTER)
	add_definitions(-DFILTER)
endif()

# Include Directories (for all targets)
include_directories(include ${CMAKE_BINARY_DIR})

# wildcard all the sources in src
file(GLOB SOURCES src/*.c src/*.cpp src/*.f95)

add_library(micropp ${SOURCES})

enable_testing()
add_subdirectory(test)
