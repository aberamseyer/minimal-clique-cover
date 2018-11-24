CC = g++
FLAGS = -std=c++11 -I boost_1_67_0

all: main

main: main.o
	${CC} ${FLAGS} main.o -o out

main.o:
	${CC} ${FLAGS} -c main.cpp

run: clean main
	./out 30.txt

mpi:
	mpicc -g -Wall ${FLAGS} -o out main.cpp

runmpi:
	mpiexec -n 4 ./out

clean:
	rm *.o
	rm out
