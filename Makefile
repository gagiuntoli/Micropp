
CC=g++
FC=gfortran
CFLAGS= -c -g -std=c++11
FFLAGS= -c -g 

#LFLAGS= -L/apps/BOOST/1.67.0/lib -lboost_program_options
#LFLAGS= -L/apps/BOOST/1.64.0_py3/INTEL/IMPI/lib -lboost_program_options
LFLAGS= -lboost_program_options

#INC= -I/apps/BOOST/1.64.0_py3/INTEL/IMPI/include
#INC=/apps/BOOST/1.67.0/include
#INC=

all: build test_1 test_2 test_3

micropp: build/main.o build/assembly.o build/micro.o build/ell.o
	$(CC) $(LFLAGS) $^ -o $@

test_1: build/test_1.o build/libmicropp.a
	$(CC) $< -o $@ $(LFLAGS) -L build -lmicropp

test_2: build/test_2.o build/libmicropp.a
	$(CC) $< -o $@ $(LFLAGS) -L build -lmicropp

test_3: build/test_3.o build/libmicropp.a
	$(FC) $< -o $@ -L build -lmicropp -lstdc++ 

build/libmicropp.a: build/assembly.o build/solve.o build/output.o  build/micro.o build/ell.o build/loc_hom.o build/wrapper.o  
	ar rcs $@ $^
    
build/main.o: src/main.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/micro.o: src/micro.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/assembly.o: src/assembly.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/ell.o: src/ell.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/test_1.o: test/test_1.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc $(INC)  -o $@

build/test_2.o: test/test_2.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc $(INC)  -o $@

build/solve.o: src/solve.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/output.o: src/output.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/loc_hom.o: src/loc_hom.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/wrapper.o: src/wrapper.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/test_3.o: test/test_3.f90
	$(FC) $(FFLAGS) $< -o $@ 

build:
	mkdir -p build

clean:
	rm -f micropp test_* build/*
