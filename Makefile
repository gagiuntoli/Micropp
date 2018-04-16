CFLAGS= -std=c++11 
LFLAGS= -lboost_program_options

micropp: build/main.o build/csr.o build/micro.o
	g++ $(CFLAGS) $(LFLAGS) build/main.o -I . -o micropp
    
build/main.o: src/main.cpp
	g++ $(CFLAGS) -c $< -I . -o build/main.o

build/csr.o: src/csr.cpp
	g++ $(CFLAGS) -c $< -I . -o build/csr.o

build/micro.o: src/micro.cpp
	g++ $(CFLAGS) -c $< -I . -o build/micro.o

clean:
	rm -f micropp build/*
