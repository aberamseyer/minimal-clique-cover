all:
	mpic++ -std=c++11 -g -O3 -lm -fopenmp graph/Graph.cpp main.cpp -o out

run: all
	mpiexec -n 24 ./out test_data/25.txt 1

	