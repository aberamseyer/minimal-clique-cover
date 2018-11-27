CC = g++
FLAGS = -std=c++11

all: main

main: main.o Graph.o
	${CC} ${FLAGS} main.o Graph.o -o out

main.o: Graph.o
	${CC} ${FLAGS} -c main.cpp

Graph.o:
	${CC} ${FLAGS} -c Graph.cpp

run: clean main
	./out test_data/13.txt

mpi:
	mpicc -g -Wall ${FLAGS} Graph.cpp main.cpp -o out

runmpi:
	mpiexec -n 4 ./out

clean:
	rm *.o
	rm out
