
cmake_minimum_required(VERSION 3.0)

project(MicroPP VERSION 0.1 LANGUAGES C CXX Fortran)
set (CMAKE_CXX_STANDARD 11)

option(OPENMP "Enables OpenMP paralelization" OFF)
if(OPENMP)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

option(TIMER "Enable time benchmarks instrumentation" OFF)
if(TIMER)
	add_definitions(-DTIMER)
endif()

option(FILTER "Enable Filtering on Average Stress and Ctan" OFF)
if(FILTER)
	add_definitions(-DFILTER)
endif()

set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")

# Include Directories (for all targets)
include_directories(include ${CMAKE_BINARY_DIR})

# wildcard all the sources in src
file(GLOB SOURCES src/*.c src/*.cpp src/*.f95)

add_library(micropp ${SOURCES})

enable_testing()
add_subdirectory(test)
