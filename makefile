CC = g++
FLAGS = -std=c++11 -I boost_1_67_0

all: main

main: main.o
	${CC} ${FLAGS} main.o -o out

main.o:
	${CC} ${FLAGS} -c main.cpp

clean:
	rm *.o
	rm out

run: clean main
	./out facebook_combined.txt

min: clean main
	./out min.txt