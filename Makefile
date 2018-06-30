#
# MicroPP : 
# Finite element library to solve microstructural problems for composite materials.
#
# Copyright (C) - 2018 - Guido Giuntoli <gagiuntoli@gmail.com>
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

CC=g++
FC=gfortran
CFLAGS= -c -std=c++11
FFLAGS= -c

ifeq ($(OPT),1)
 FFLAGS += -O3
 CFLAGS += -O3
else
 FFLAGS += -g 
 CFLAGS += -g
endif

all: build test_3 test_7 test_8

lib: build/libmicropp.a

test_3: build/test_3.o build/libmicropp.a
	$(FC) $< -o $@ -L build -lmicropp -lstdc++ 

test_7: build/test_7.o build/libmicropp.a
	$(CC) $< -o $@ -L build -lmicropp 

test_8: build/test_8.o build/libmicropp.a
	$(CC) $< -o $@ -L build -lmicropp 

build/libmicropp.a: build/assembly.o build/solve.o build/output.o  build/micro.o build/ell.o build/loc_hom.o build/wrapper.o  
	ar rcs $@ $^
    
build/%.o: test/%.f90
	$(FC) $(FFLAGS) $< -o $@ 

build/%.o: test/%.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc $(INC)  -o $@

build/%.o: src/%.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc $(INC)  -o $@


.PHONY: clean cleanouts

clean:
	rm -f micropp test_* build/*

cleanouts:
	rm -f micropp_* *.dat

build:
	mkdir -p build
