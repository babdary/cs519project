compile:
	g++ generation_code.c -o generate_graph 
	g++ -std=c++11 -fopenmp -o floyd_warshall FloydWarshall_OMP.cc

run:
	./generate_graph 6000 12000
	./floyd_warshall inputMatrix.txt 12
