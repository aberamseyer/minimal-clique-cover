CC = g++
FLAGS = -std=c++11 -g -Wall -O3

all: main

main: main.o Graph.o
	${CC} ${FLAGS} main.o Graph.o -o out

main.o: Graph.o
	${CC} ${FLAGS} -c main.cpp

Graph.o:
	${CC} ${FLAGS} -c Graph.cpp

run: clean main
	clear
	./out test_data/13.txt

openmp:
	${CC} ${FLAGS} -fopenmp -pg Graph.cpp main.cpp -o out

mpi:
	mpic++ ${FLAGS} -lm Graph.cpp main.cpp -o out

runmpi:
	mpiexec -n 4 ./out

clean:
	rm *.o
	rm out
