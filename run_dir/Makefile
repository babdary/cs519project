#batch_mpi.sh  fw_par_block.c  fw_par_row  fw_par_row_block_cyclic  fw_par_row_block_cyclic.c  fw_par_row.c  fw_ser.c  generate_graph  generation_code.c  inputMatrix.txt  Makefile

compile:
	g++ generation_code.c -o generate_graph 
	mpicxx -o fw_par_row fw_par_row.c
	mpicxx -o fw_par_row_block_cyclic fw_par_row_block_cyclic.c
	mpicxx -o fw_par_block fw_par_block.c
	g++ -o fw_ser fw_ser.c

run:
	./generate_graph 4 8
