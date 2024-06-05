/*ROW parallelization NOT VERY OPTIMIZED THO xd*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void floyd_warshall(int num_of_vertices, int local_iter, double *local_adj_matrix, int my_rank, int comm_sz){
    
    int i,j,k;
    double dist;
    double *row_k =NULL;

    for(k = 0; k < num_of_vertices; k++){

        int root = k/local_iter;
        row_k = (double *)malloc(num_of_vertices*sizeof(double));

        if(my_rank == root){
            for(int j = 0; j < num_of_vertices; j++){
                row_k[j] = local_adj_matrix[(k % local_iter)* num_of_vertices +j];
            }
        }
        MPI_Bcast(row_k, num_of_vertices,MPI_DOUBLE, root,MPI_COMM_WORLD);
        for(i = 0; i < local_iter; i++){
            for(j = 0; j < num_of_vertices; j++){
                dist = local_adj_matrix[i * num_of_vertices + k] + row_k[j];
                if(local_adj_matrix[i * num_of_vertices + j] > dist){
                    local_adj_matrix[i * num_of_vertices + j] = dist;
                }
            }
        }
    }
    free(row_k);
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

    int num_of_vertices;
    char *input_file = argv[1];
    double *matrix = NULL;

    if (my_rank == 0)
    {
        if (argc < 2)
        {
            printf("Provide ./executable <input file name>");
        }
        file = fopen(input_file, "r");
        if (file == NULL)
        {
            printf("Cannot open input file!!!");
        }
        fscanf(file, "%d", &num_of_vertices);
        //fprintf(stderr, "number of vertices = %d\n", num_of_vertices);

 
        matrix = (double *)malloc(num_of_vertices * num_of_vertices * sizeof(double));
        for (int ii = 0; ii < num_of_vertices; ii++)
        {
            for (int jj = 0; jj < num_of_vertices; jj++)
            {
                fscanf(file, "%lf", &matrix[ii * num_of_vertices + jj]);
            }
        }
        fclose(file);
        // for (int i = 0; i < num_of_vertices; i++)
        // {
        //     for (int j = 0; j < num_of_vertices; j++)
        //     {
        //         printf("row %d = %lf\n", i, matrix[i * num_of_vertices + j]);
        //     }
        // }
    }

    
    MPI_Bcast(&num_of_vertices, 1, MPI_INT, 0, MPI_COMM_WORLD); // let each proc know the num of vertices
    int local_iter = num_of_vertices / comm_sz;                 // number of rows each proc will handle
    // printf("local_iter: %d\n", local_iter);
    /*Read Adjacency Matrix of the Weighted Directed Graph from Input File*/
    double *local_adj_matrix = (double *)malloc(local_iter * num_of_vertices * sizeof(double));
  
    MPI_Scatter(matrix, local_iter * num_of_vertices, MPI_DOUBLE, local_adj_matrix, local_iter * num_of_vertices, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    
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
                printf("| %lf ", matrix[i* num_of_vertices +j]);
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