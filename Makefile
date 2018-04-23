
CFLAGS= -c -std=c++11 -O3
#LFLAGS= -L/apps/BOOST/1.67.0/lib -lboost_program_options
LFLAGS= -lboost_program_options
#INC=/apps/BOOST/1.67.0/include
INC=

micropp: build/main.o build/csr.o build/micro.o
	g++ $(LFLAGS) $^ -o micropp
    
build/main.o: src/main.cpp
	g++ $(CFLAGS) $< -I ./inc -o build/main.o

build/micro.o: src/micro.cpp
	g++ $(CFLAGS) $< -I ./inc -o build/micro.o

build/assembly.o: src/assembly.cpp
	g++ $(CFLAGS) $< -I ./inc -o build/assembly.o

build/csr.o: src/csr.cpp
	g++ $(CFLAGS) $< -I ./inc -o build/csr.o

build/ell.o: src/ell.cpp
	g++ $(CFLAGS) $< -I ./inc -o build/micro.o

clean:
	rm -f micropp build/*
