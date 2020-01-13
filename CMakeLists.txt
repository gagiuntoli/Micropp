
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
	file(GLOB SOURCES src/*.c src/*.cpp src/*.f95 src/cuda/*.cu)
        set_source_files_properties(src/common.cpp PROPERTIES LANGUAGE CUDA)
elseif(ENABLE_OPENACC)
	file(GLOB SOURCES 
		src/material.c 
		src/material.cpp
		src/material.f95 
		src/micropp.cpp 
		src/micropp_c.cpp 
		src/micropp.f95
		src/openacc/assembly.cpp 
		src/openacc/ell.cpp 
		src/ell-common.cpp 
		src/instrument.cpp 
		src/homogenize.cpp 
		src/common.cpp 
		src/solve.cpp 
		src/micro3D.cpp 
		src/output.cpp 
		)
else()
	file(GLOB SOURCES 
		src/material.c 
		src/material.cpp
		src/material.f95 
		src/micropp.cpp 
		src/micropp_c.cpp 
		src/micropp.f95
		src/assembly.cpp 
		src/ell.cpp 
		src/ell-common.cpp 
		src/instrument.cpp 
		src/homogenize.cpp 
		src/common.cpp 
		src/solve.cpp 
		src/micro3D.cpp 
		src/output.cpp 
		)
endif()

add_library(micropp ${SOURCES})

if(ENABLE_CUDA)
set_property(TARGET micropp PROPERTY CUDA_SEPARABLE_COMPILATION ON)
endif()

enable_testing()
add_subdirectory(test)
