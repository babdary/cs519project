/*ROW parallelization NOT VERY OPTIMIZED THO xd*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



void floyd warshall(int num_of_vertices, int local_iter, double *local_adj_matrix, int my_rank, int comm_sz){
    
    for(int n = 0; n < num_of_vertices; n++){

        int root = n/local_iter;
        double *node_path = (double *)malloc(num_of_vertices*sizeof(double));

        if(my_rank = root){
            memcpy(node_path, &local_adj_matrix[(n%local_iter)*num_of_vertices], num_of_vertices*sizeof(double));
        }
    }

}

int main(int argc, char *argv[])
{

    int my_rank, comm_sz;
    FILE *file;

    /*MPI initialization*/
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    int num_of_vertices;
    char *input_file = argv[1];

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
        fprintf(stderr, "number of vertices = %d\n", num_of_vertices);

        double *matrix = (double *)malloc(num_of_vertices * num_of_vertices * sizeof(double));
        for (int ii = 0; ii < num_of_vertices; ii++)
        {
            for (int jj = 0; jj < num_of_vertices; jj++)
            {
                fscanf(file, "%lf", &matrix[ii * num_of_vertices + jj]);
            }
        }
        fclose(file);
        for (int i = 0; i < num_of_vertices; i++)
        {
            for (int j = 0; j < num_of_vertices; j++)
            {
                printf("row %d = %lf\n", i, matrix[i * num_of_vertices + j]);
            }
        }
    }

    
    MPI_Bcast(&num_of_vertices, 1, MPI_INT, 0, MPI_COMM_WORLD); // let each proc know the num of vertices
    int local_iter = num_of_vertices / comm_sz;                 // number of rows each proc will handle
    printf("local_iter: %d\n", local_iter);
    /*Read Adjacency Matrix of the Weighted Directed Graph from Input File*/
    double *local_adj_matrix = (double *)malloc(local_iter * num_of_vertices * sizeof(double));
   
    MPI_Scatter(matrix, local_iter * number_of_vertices, MPI_DOUBLE, local_adj_matrix, local_iter * number_of_vertices, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;

}