/*ROW parallelization NOT VERY OPTIMIZED THO xd*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define min(x, y) (((x) < (y)) ? (x) : (y))


#define INF 1e10

void floyd_warshall(int num_of_vertices, int local_iter, double *local_adj_matrix, int my_rank, int comm_sz){
    
    int i,j,k;
    double dist;
    double *k_row = (double *)malloc(num_of_vertices*sizeof(double));

    for(k = 0; k < num_of_vertices; k++){

        int root = k/local_iter;

        if(my_rank == root){
            for(int j = 0; j < num_of_vertices; j++){
                k_row[j] = local_adj_matrix[(k % local_iter)* num_of_vertices +j];
            }
        }
        MPI_Bcast(k_row, num_of_vertices,MPI_DOUBLE, root, MPI_COMM_WORLD);
        for(i = 0; i < local_iter; i++){
            for(j = 0; j < num_of_vertices; j++){
                dist = local_adj_matrix[i * num_of_vertices + k] + k_row[j];
                if(local_adj_matrix[i * num_of_vertices + j] > dist){
                    local_adj_matrix[i * num_of_vertices + j] = dist;
                }
            }
        }
    }
    free(k_row);
}

int main(int argc, char *argv[])
    {

        clock_t fw_start, fw_end, start, end;

        double wall_clock_time, fw_time;

        start = clock();

        int my_rank, comm_sz;
        FILE *file;

        /*MPI initialization*/
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

        int num_of_vertices, number_of_edges;
        char *input_file = argv[1];
        double *matrix = NULL;

        if (my_rank == 0) {
        file = fopen(input_file, "r");
        if (file == NULL) {
            printf("Cannot open input file!!!\n");
            MPI_Finalize();
            return -1;
        }

        fscanf(file, "%d %d", &num_of_vertices, &number_of_edges);

        // Allocate global data
        matrix = (double *)malloc(num_of_vertices * num_of_vertices * sizeof(double));
        for (int i = 0; i < num_of_vertices; i++) {
            for (int j = 0; j < num_of_vertices; j++) {
                matrix[i*num_of_vertices + j] = (i == j) ? 0 : INF;
            }
        }

        int source, target;
        double weight;
        // Store the weights in a matrix
        for (int i = 0; i < number_of_edges; i++) {
            fscanf(file, "%d %d %lf", &source, &target, &weight);
            source--; target--; // Convert 1-based to 0-based index
            matrix[source*num_of_vertices + target ] = min(matrix[source*num_of_vertices + target], weight);
        }
        fclose(file);
    }

    MPI_Bcast(&num_of_vertices, 1, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast the number of vertices

    // Calculate local_iter and allocate memory for local_adj_matrix
    int local_iter = num_of_vertices / comm_sz; // Number of rows each process will handle
    double *local_adj_matrix = (double *)malloc(local_iter * num_of_vertices * sizeof(double));

    // Scatter the matrix data from process 0 to all other processes
    MPI_Scatter(matrix, local_iter * num_of_vertices, MPI_DOUBLE, local_adj_matrix, local_iter * num_of_vertices, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Now you can perform Floyd Warshall algorithm on each process with the received local_adj_matrix

    

    /*Call Floyd Warshall*/
    fw_start = clock();
    floyd_warshall(num_of_vertices,local_iter,local_adj_matrix,my_rank,comm_sz);
    fw_end = clock();
    fw_time = ((double)(fw_end-fw_start)) /CLOCKS_PER_SEC;



    /*Gather results*/
    MPI_Gather(local_adj_matrix, local_iter * num_of_vertices, MPI_DOUBLE, matrix, local_iter * num_of_vertices, MPI_DOUBLE, 0 , MPI_COMM_WORLD);


    if(my_rank==0){
        int i,j;
        for(i = 0; i < num_of_vertices; i++){
            for(j = 0; j < num_of_vertices; j++){
                if(matrix[i* num_of_vertices +j]==INF){
                    printf("%6s", "INF");
                }
                    printf("| %6.2f ", matrix[i* num_of_vertices +j]);
            }
            printf("|\n");
        }
        free(matrix);
    }
    free(local_adj_matrix);
    MPI_Finalize();


    end = clock();
    wall_clock_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time Elapsed:%f\n", wall_clock_time);
    printf("Floyd Warshall Algorithm Elapsed: %f\n", fw_time);
    printf("------------------------------------------------------\n");
    return 0;

}