#
#  This source code is part of MicroPP: a finite element library
#  to solve microstructural problems for composite materials.
#
#  Copyright (C) - 2018 - Guido Giuntoli <gagiuntoli@gmail.com>
#                         Jimmy Aguilar Mena <kratsbinovish@gmail.com>
#                         Judicaël Grasset <judicael.grasset@stfc.ac.uk>
#                         Alejandro Figueroa <afiguer7@maisonlive.gmu.edu>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published
#  by the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

cmake_minimum_required(VERSION 3.20)

project(
  MicroPP
  VERSION 0.1
  LANGUAGES C CXX Fortran)
set(CMAKE_CXX_STANDARD 11)

set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")

option(ENABLE_CUDA "Enable CUDA" OFF)
if(ENABLE_CUDA)
	enable_language(CUDA)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_CUDA")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -D_CUDA")
endif()

option(ENABLE_OPENMP "Enable OpenMP" OFF)
if(OPENMP)
  find_package(OpenMP REQUIRED)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

option(ENABLE_OPENACC "Enable OpenACC" OFF)
if(ENABLE_OPENACC)
  find_package(OpenACC REQUIRED)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenACC_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenACC_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenACC_Fortran_FLAGS}")
endif()

option(ENABLE_TIMER "Enable instrumentation" OFF)
if(ENABLE_TIMER)
  add_definitions(-DTIMER)
endif()

# Include Directories (for all targets)
include_directories(include ${CMAKE_BINARY_DIR})

file(GLOB SOURCES_COMMON_CPP
     src/material.cpp
     src/micropp.cpp
     src/micropp_c.cpp 
     src/average.cpp 
     src/update.cpp 
     src/ell-common.cpp 
     src/homogenize.cpp 
     src/common.cpp 
     src/solve.cpp 
     src/micro3D.cpp 
     src/output.cpp 
)

# wildcard all the sources in src
if(ENABLE_CUDA)
	file(GLOB SOURCES 
	        ${SOURCES_COMMON_CPP}
		src/material.c 
		src/material.f95 
		src/micropp.f95
		src/cuda/assembly.cu 
		src/cuda/ell.cu 
		src/cuda/init.cu 
		src/instrument.cpp 
		)
        set_source_files_properties(${SOURCES_COMMON_CPP} PROPERTIES LANGUAGE CUDA)
        set_source_files_properties(src/material.c PROPERTIES LANGUAGE CUDA)
elseif(ENABLE_OPENACC)
	file(GLOB SOURCES 
	        ${SOURCES_COMMON_CPP}
		src/material.c 
		src/material.f95 
		src/micropp.f95
		src/openacc/assembly.cpp 
		src/openacc/ell.cpp 
		src/instrument.cpp 
		)
else()
	file(GLOB SOURCES 
	        ${SOURCES_COMMON_CPP}
		src/material.c 
		src/material.f95 
		src/micropp.f95
		src/assembly.cpp 
		src/average.cpp 
		src/update.cpp 
		src/ell.cpp 
		src/instrument.cpp 
		)
endif()

add_library(micropp ${SOURCES})

if(ENABLE_CUDA)
  set_property(TARGET micropp PROPERTY CUDA_SEPARABLE_COMPILATION ON)
endif()

enable_testing()
add_subdirectory(test)
