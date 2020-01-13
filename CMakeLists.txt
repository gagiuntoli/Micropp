
cmake_minimum_required(VERSION 3.0)

project(MicroPP VERSION 0.1 LANGUAGES C CXX Fortran)
set (CMAKE_CXX_STANDARD 11)

set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")

option(ENABLE_CUDA "Enable CUDA" OFF)
if(ENABLE_CUDA)
	enable_language(CUDA)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_CUDA")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -DCUDA")
endif()

option(ENABLE_OPENMP "Enables OpenMP paralelization" OFF)
if(ENABLE_OPENMP)
	if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
		set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp ")
	elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mp ")
		set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mp ")
	endif()
endif()

option(ENABLE_OPENACC "Enable OpenACC" OFF)
if(ENABLE_OPENACC)
	#add_compile_options(-acc -Minfo=accel)
	if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenacc")
		set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenacc")
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenacc")
	elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc -Minfo=accel")
		set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -acc -Minfo=accel")
		if(CMAKE_BUILD_TYPE EQUAL "Debug" OR CMAKE_BUILD_TYPE EQUAL "RelWithDebInfo")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Mprof=ccff")
		endif()
	endif()
endif()


option(ENABLE_TIMER "Enable time benchmarks instrumentation" OFF)
if(ENABLE_TIMER)
	add_definitions(-DTIMER)
endif()

# Include Directories (for all targets)
include_directories(include ${CMAKE_BINARY_DIR})

# wildcard all the sources in src
if(ENABLE_CUDA)
	file(GLOB SOURCES src/*.c src/*.cpp src/*.f95 src/*.cu)
        set_source_files_properties(src/common.cpp PROPERTIES LANGUAGE CUDA)
else()
	file(GLOB SOURCES src/*.c src/*.cpp src/*.f95)
endif()

add_library(micropp ${SOURCES})

if(ENABLE_CUDA)
set_property(TARGET micropp PROPERTY CUDA_SEPARABLE_COMPILATION ON)
endif()

enable_testing()
add_subdirectory(test)
