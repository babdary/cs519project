compile:
	g++ generation_code.c -o generate_graph 
	mpicxx -o floyd_warshall_mpi FloydWarshall_MPI.cc

run:
	./generate_graph 4 8
	mpirun -np 2 ./floyd_warshall_mpi  inputMatrix.txt 
