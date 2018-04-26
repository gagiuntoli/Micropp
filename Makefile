
CFLAGS= -c -g -std=c++11 -O3
#LFLAGS= -L/apps/BOOST/1.67.0/lib -lboost_program_options
LFLAGS= -lboost_program_options
#INC=/apps/BOOST/1.67.0/include
INC=

all: test_1

micropp: build/main.o build/assembly.o build/csr.o build/micro.o build/ell.o
	g++ $(LFLAGS) $^ -o $@

test_1: build/test_1.o build/libmicropp.a
	g++ $< -o $@ $(LFLAGS) -L build -lmicropp

build/libmicropp.a: build/assembly.o build/solve.o build/csr.o build/micro.o build/ell.o 
	ar rcs $@ $^
    
build/main.o: src/main.cpp inc/micro.h
	g++ $(CFLAGS) $< -I ./inc -o $@

build/micro.o: src/micro.cpp inc/micro.h
	g++ $(CFLAGS) $< -I ./inc -o $@

build/assembly.o: src/assembly.cpp inc/micro.h
	g++ $(CFLAGS) $< -I ./inc -o $@

build/csr.o: src/csr.cpp inc/micro.h
	g++ $(CFLAGS) $< -I ./inc -o $@

build/ell.o: src/ell.cpp inc/micro.h
	g++ $(CFLAGS) $< -I ./inc -o $@

build/test_1.o: src/test_1.cpp inc/micro.h
	g++ $(CFLAGS) $< -I ./inc -o $@

build/solve.o: src/solve.cpp inc/micro.h
	g++ $(CFLAGS) $< -I ./inc -o $@

clean:
	rm -f micropp build/*
