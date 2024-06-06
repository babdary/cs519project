compile:
	g++ generation_code.c -o generate_graph 
	mpicxx -o fw_par_row fw_par_row.c
	mpicxx -o fw_par_row_block_cyclic fw_par_row_block_cyclic.c
	mpicxx -o fw_par_block fw_par_block.c
	g++ -o fw_ser fw_ser.c
	mpicxx -o fw_par_3 FWPar3.c

generate_inputs:
	./generate_graph 40 80
