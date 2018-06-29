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

#LFLAGS= -L/apps/BOOST/1.67.0/lib -lboost_program_options
#LFLAGS= -L/apps/BOOST/1.64.0_py3/INTEL/IMPI/lib -lboost_program_options
LFLAGS += -lboost_program_options

#INC= -I/apps/BOOST/1.64.0_py3/INTEL/IMPI/include
#INC=/apps/BOOST/1.67.0/include
#INC=

all: build test_1 test_3 test_5 test_6 test_7 test_8

lib: build/libmicropp.a

test_1: build/test_1.o build/libmicropp.a
	$(CC) $< -o $@ $(LFLAGS) -L build -lmicropp

test_2: build/test_2.o build/libmicropp.a
	$(CC) $< -o $@ $(LFLAGS) -L build -lmicropp

test_3: build/test_3.o build/libmicropp.a
	$(FC) $< -o $@ -L build -lmicropp -lstdc++ 

test_4: build/test_4.o build/libmicropp.a
	$(CC) $< -o $@ $(LFLAGS) -L build -lmicropp 

test_5: build/test_5.o build/libmicropp.a
	$(CC) $< -o $@ $(LFLAGS) -L build -lmicropp 

test_6: build/test_6.o build/libmicropp.a
	$(CC) $< -o $@ $(LFLAGS) -L build -lmicropp 

test_7: build/test_7.o build/libmicropp.a
	$(CC) $< -o $@ -L build -lmicropp 

test_8: build/test_8.o build/libmicropp.a
	$(CC) $< -o $@ -L build -lmicropp 

build/libmicropp.a: build/assembly.o build/solve.o build/output.o  build/micro.o build/ell.o build/loc_hom.o build/wrapper.o  
	ar rcs $@ $^
    
build/main.o: src/main.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/micro.o: src/micro.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/assembly.o: src/assembly.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/ell.o: src/ell.cpp inc/micro.h inc/ell.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/test_1.o: test/test_1.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc $(INC)  -o $@

build/test_2.o: test/test_2.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc $(INC)  -o $@

build/test_3.o: test/test_3.f90
	$(FC) $(FFLAGS) $< -o $@ 

build/test_4.o: test/test_4.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc $(INC)  -o $@

build/test_5.o: test/test_5.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc $(INC)  -o $@

build/test_6.o: test/test_6.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc $(INC)  -o $@

build/test_7.o: test/test_7.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc $(INC)  -o $@

build/test_8.o: test/test_8.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc $(INC)  -o $@

build/solve.o: src/solve.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/output.o: src/output.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/loc_hom.o: src/loc_hom.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/wrapper.o: src/wrapper.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build:
	mkdir -p build

clean:
	rm -f micropp test_* build/*
