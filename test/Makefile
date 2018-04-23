
OPT = 0
ifeq ($(OPT),1)
CFLAGS=-O3
else
CFLAGS=-g -c -O0 -std=c++11
endif

all: ell_solvers

ell_solvers: build/ell_solvers.o build/ell.o
	g++ -o $@ $^ 

build/ell_solvers.o: ell_solvers.cpp
	g++ ${CFLAGS} -o $@ $< -I ../inc

build/ell.o: ../src/ell.cpp ../inc/ell.h
	g++ ${CFLAGS} -o $@ $< -I ../inc

clean:	    
	rm -f build/* *.o ell_solvers.o

run:
	./ell_solvers
