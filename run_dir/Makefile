compile:
	g++ generation_code.c -o generate_graph 
	mpicxx -o fw_par_row fw_par_row.c
	mpicxx -o fw_par_row_block_cyclic fw_par_row_block_cyclic.c
	mpicxx -o fw_par_block fw_par_block.c
	g++ -o fw_ser fw_ser.c
	mpicxx -o fw_par_3 FWPar3.c

generate_inputs:
	./generate_graph 10000 20000

run_serial:
	./fw_ser inputMatrix.txt

run_row:
	mpirun -np 4 ./fw_par_row inputMatrix.txt

run_row_cyclic:
	mpirun -np 4 ./fw_par_row_block_cyclic inputMatrix.txt

run_par_3:
	mpirun -np 4 ./fw_par_3 inputMatrix.txt

run_par_block:
	mpirun -np 4 ./fw_par_block inputMatrix.txt
