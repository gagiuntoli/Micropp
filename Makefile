CFLAGS= -c -g -std=c++11 
LFLAGS= -lboost_program_options

micropp: build/main.o build/csr.o build/micro.o
	g++ $(LFLAGS) $^ -o micropp
    
build/main.o: src/main.cpp
	g++ $(CFLAGS) $< -I inc/ -o build/main.o

build/csr.o: src/csr.cpp
	g++ $(CFLAGS) $< -I inc/ -o build/csr.o

build/micro.o: src/micro.cpp
	g++ $(CFLAGS) $< -I inc/ -o build/micro.o

clean:
	rm -f micropp build/*
