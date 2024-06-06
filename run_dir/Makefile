compile:
	g++ generation_code.c -o generate_graph 
	mpicxx -o fw_par_row fw_par_row.c
	mpicxx -o fw_par_row_block_cyclic fw_par_row_block_cyclic.c

run:
	./generate_graph 4 8
