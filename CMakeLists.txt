cmake_minimum_required(VERSION 3.0)

project(
  MicroPP
  VERSION 0.1
  LANGUAGES C CXX Fortran)
set(CMAKE_CXX_STANDARD 11)

set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")

option(OPENMP "Enables OpenMP paralelization" OFF)
if(OPENMP)
  find_package(OpenMP REQUIRED)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

option(OPENACC "Enable OpenACC" OFF)
if(OPENACC)
  find_package(OpenACC REQUIRED)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenACC_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenACC_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenACC_Fortran_FLAGS}")
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    if(CMAKE_BUILD_TYPE EQUAL "Debug" OR CMAKE_BUILD_TYPE EQUAL "RelWithDebInfo")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Mprof=ccff")
    endif()
  endif()
endif()

option(TIMER "Enable time benchmarks instrumentation" OFF)
if(TIMER)
  add_definitions(-DTIMER)
endif()

# Include Directories (for all targets)
include_directories(include ${CMAKE_BINARY_DIR})

# wildcard all the sources in src
file(GLOB SOURCES src/*.c src/*.cpp src/*.f95)

add_library(micropp ${SOURCES})

enable_testing()
add_subdirectory(test)
