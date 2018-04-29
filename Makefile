CC=g++
CFLAGS= -c -g -std=c++11
#LFLAGS= -L/apps/BOOST/1.67.0/lib -lboost_program_options
LFLAGS= -lboost_program_options
#INC=/apps/BOOST/1.67.0/include
INC=

all: test_1 test_2

micropp: build/main.o build/assembly.o build/micro.o build/ell.o
	$(CC) $(LFLAGS) $^ -o $@

test_1: build/test_1.o build/libmicropp.a
	$(CC) $< -o $@ $(LFLAGS) -L build -lmicropp

test_2: build/test_2.o build/libmicropp.a
	$(CC) $< -o $@ $(LFLAGS) -L build -lmicropp

build/libmicropp.a: build/assembly.o build/solve.o build/output.o  build/micro.o build/ell.o build/loc_hom.o 
	ar rcs $@ $^
    
build/main.o: src/main.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/micro.o: src/micro.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/assembly.o: src/assembly.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/ell.o: src/ell.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/test_1.o: src/test_1.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/test_2.o: src/test_2.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/solve.o: src/solve.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/output.o: src/output.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

build/loc_hom.o: src/loc_hom.cpp inc/micro.h
	$(CC) $(CFLAGS) $< -I ./inc -o $@

clean:
	rm -f micropp test_1 build/*
