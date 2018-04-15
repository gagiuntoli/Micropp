micropp: build/main.o build/csr.o build/micro.o
	g++ build/main.o -I . -o micropp
    
build/main.o: src/main.cpp
	g++ -c $< -I . -o build/main.o

build/csr.o: src/csr.cpp
	g++ -c $< -I . -o build/csr.o

build/micro.o: src/micro.cpp
	g++ -c $< -I . -o build/micro.o

clean:
	rm -f micropp build/*
