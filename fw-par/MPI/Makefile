compile:
	g++ generation_code.c -o generate_graph 
	mpicxx -o floyd_warshall_mpi FloydWarshall_MPI.cc

run:
	./generate_graph 8 30
	mpirun -np 12 ./floyd_warshall_mpi  inputMatrix.txt output.txt 1
