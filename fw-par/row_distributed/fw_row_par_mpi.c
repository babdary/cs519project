/*ROW parallelization NOT VERY OPTIMIZED THO xd*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char *argv[])
{

    int my_rank, comm_sz;
    FILE *file;

    /*MPI initialization*/
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    int number_of_vertices;
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
        fscanf(file, "%d", &number_of_vertices);
        fprintf(stderr, "number of vertices = %d\n", number_of_vertices);

        double *matrix = (double *)malloc(number_of_vertices * number_of_vertices * sizeof(double));
        for (int ii = 0; ii < number_of_vertices; ii++)
        {
            for (int jj = 0; jj < number_of_vertices; jj++)
            {
                fscanf(file, "%lf", &matrix[ii * number_of_vertices + jj]);
            }
        }
        fclose(file);
        for (int i = 0; i < number_of_vertices; i++)
        {
            for (int j = 0; j < number_of_vertices; j++)
            {
                printf("row %d = %lf\n", i, matrix[i * number_of_vertices + j]);
            }
        }
    }

    
    MPI_Bcast(&number_of_vertices, 1, MPI_INT, 0, MPI_COMM_WORLD); // let each proc know the num of vertices
    int local_iter = number_of_vertices / comm_sz;                 // number of rows each proc will handle
    printf("local_iter: %d\n", local_iter);

    double *local_adj_matrix = (double *)malloc(local_iter * number_of_vertices * sizeof(double));
    /*Read Adjacency Matrix of the Weighted Directed Graph from Input File*/

    MPI_Finalize();
    return 0;
}