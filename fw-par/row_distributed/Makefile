compile:
	g++ generation_code.c -o generate_graph 
	mpicc  -lm fw_row_par_mpi_sparse.c -o fw_row_par_mpi_sparse

run: compile
	./generate_graph 6 12
	mpirun -np 2 ./fw_row_par_mpi_sparse  inputMatrix.txt 
