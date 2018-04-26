
CFLAGS= -c -g -std=c++11 -O3
#LFLAGS= -L/apps/BOOST/1.67.0/lib -lboost_program_options
LFLAGS= -lboost_program_options
#INC=/apps/BOOST/1.67.0/include
INC=

micropp: build/main.o build/assembly.o build/csr.o build/micro.o build/ell.o
	g++ $(LFLAGS) $^ -o $@

test_1: build/test_1.o build/libmicropp.a
	g++ $< -o $@ $(LFLAGS) -L build -lmicropp

build/libmicropp.a: build/assembly.o build/csr.o build/micro.o build/ell.o
	ar rcs $@ $^
    
build/main.o: src/main.cpp
	g++ $(CFLAGS) $< -I ./inc -o $@

build/micro.o: src/micro.cpp
	g++ $(CFLAGS) $< -I ./inc -o $@

build/assembly.o: src/assembly.cpp
	g++ $(CFLAGS) $< -I ./inc -o $@

build/csr.o: src/csr.cpp
	g++ $(CFLAGS) $< -I ./inc -o $@

build/ell.o: src/ell.cpp
	g++ $(CFLAGS) $< -I ./inc -o $@

build/test_1.o: src/test_1.cpp
	g++ $(CFLAGS) $< -I ./inc -o $@

clean:
	rm -f micropp build/*
