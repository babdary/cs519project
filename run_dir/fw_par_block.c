#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

#define INF 99999.0
#define MIN(A, B) (A < B) ? A : B

void printSolution(double *dist, int V)
{
    FILE *outfile = fopen("out_block.txt", "w");
    for (int i = 0; i < V; i++)
    {
        for (int j = 0; j < V; j++)
        {
            if (dist[i * V + j] == INF)
                fprintf(outfile, "| %6s ", "INF");
            else
                fprintf(outfile, "| %6.2f ", dist[i * V + j]);
        }
        fprintf(outfile, "\n");
    }
    fclose(outfile);
}

int main(int argc, char *argv[])
{
    struct timespec start, end, floyd_start, floyd_end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    MPI_Init(&argc, &argv);

    FILE *file;
    int V, E, my_rank, comm_sz, local_n, col_rank, row_rank, block = 0;
    int dims[2], periods[2], coords[2];
    double *graph = NULL;

    MPI_Comm comm_row, comm_col, comm_block;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    int square_length = (int)sqrt(comm_sz);

    dims[0] = dims[1] = square_length;
    periods[0] = periods[1] = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm_block);

    MPI_Cart_coords(comm_block, my_rank, 2, coords);

    coords[0] = 0;
    coords[1] = 1;
    MPI_Cart_sub(comm_block, coords, &comm_row);
    coords[0] = 1;
    coords[1] = 0;
    MPI_Cart_sub(comm_block, coords, &comm_col);

    MPI_Comm_rank(comm_row, &row_rank);
    MPI_Comm_rank(comm_col, &col_rank);

    if (my_rank == 0)
    {
        file = fopen(argv[1], "r");
        fscanf(file, "%d", &V);
        fscanf(file, "%d", &E);

        if( square_length * square_length != comm_sz || V % square_length != 0)
        {
            printf("Invalid number of processors or matrix size.\n");
            MPI_Abort(MPI_COMM_WORLD, 0);
            exit(1);
        }

        graph = (double *)malloc(V * V * sizeof(double));
        for (int i = 0; i < V; i++)
        {
            for (int j = 0; j < V; j++)
            {
                if (i == j)
                    graph[i * V + j] = 0.0;
                else
                    graph[i * V + j] = INF;
            }
        }

        int u, v;
        double w;
        for (int i = 0; i < E; i++)
        {
            fscanf(file, "%d %d %lf", &u, &v, &w);
            graph[(u - 1) * V + (v - 1)] = w;
        }

        fclose(file);
    }

    MPI_Bcast(&V, 1, MPI_INT, 0, MPI_COMM_WORLD);

    local_n = V / square_length;
    clock_gettime(CLOCK_MONOTONIC, &floyd_start);

    double *sendbuf = NULL;
    int *sendcounts = NULL;
    int *displs = NULL;

    if (my_rank == 0)
    {
        sendbuf = (double *)malloc(V * V * sizeof(double));
        sendcounts = (int *)malloc(comm_sz * sizeof(int));
        displs = (int *)malloc(comm_sz * sizeof(int));

        for (int i = 0; i < comm_sz; i++)
        {
            sendcounts[i] = local_n * local_n;
            displs[i] = i * local_n * local_n;
        }

        int index = 0;
        for (int k = 0; k < square_length; k++)
        {
            for (int l = 0; l < square_length; l++)
            {
                for (int i = 0; i < local_n; i++)
                {
                    for (int j = 0; j < local_n; j++)
                    {
                        sendbuf[index++] = graph[(k * local_n + i) * V + (l * local_n + j)];
                    }
                }
            }
        }
    }

    double *recvbuf = (double *)malloc(local_n * local_n * sizeof(double));

    MPI_Scatterv(sendbuf, sendcounts, displs, MPI_DOUBLE, recvbuf, local_n * local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double *local_col = (double *)malloc(local_n * sizeof(double));
    double *local_row = (double *)malloc(local_n * sizeof(double));

    for (int k = 0; k < V; k++)
    {
        int index = k % local_n;

        if (k > 0 && index == 0)
        {
            block++;
        }

        if (col_rank == block)
        {
            for (int i = 0; i < local_n; i++)
            {
                local_row[i] = recvbuf[index * local_n + i];
            }
        }

        if (row_rank == block)
        {
            for (int i = 0; i < local_n; i++)
            {
                local_col[i] = recvbuf[i * local_n + index];
            }
        }

        MPI_Bcast(local_row, local_n, MPI_DOUBLE, block, comm_col);
        MPI_Bcast(local_col, local_n, MPI_DOUBLE, block, comm_row);

        for (int i = 0; i < local_n; i++)
        {
            for (int j = 0; j < local_n; j++)
            {
                recvbuf[i * local_n + j] = MIN(recvbuf[i * local_n + j], local_col[i] + local_row[j]);
            }
        }
    }

    MPI_Gatherv(recvbuf, local_n * local_n, MPI_DOUBLE, graph, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        clock_gettime(CLOCK_MONOTONIC, &floyd_end);
        double *result = (double *)malloc(V * V * sizeof(double));

        int index = 0;
        for (int k = 0; k < square_length; k++)
        {
            for (int l = 0; l < square_length; l++)
            {
                for (int i = 0; i < local_n; i++)
                {
                    for (int j = 0; j < local_n; j++)
                    {
                        result[(k * local_n + i) * V + (l * local_n + j)] = graph[index++];
                    }
                }
            }
        }

      
        printSolution(result, V);

        clock_gettime(CLOCK_MONOTONIC, &end);

        double total_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
        double floyd_time = (floyd_end.tv_sec - floyd_start.tv_sec) + (floyd_end.tv_nsec - floyd_start.tv_nsec) / 1e9;

        printf("Total Time measured: %.3f seconds.\n", total_time);
        printf("Floyd Time measured: %.3f seconds.\n", floyd_time);

        free(result);
        free(graph);
        free(sendbuf);
        free(sendcounts);
        free(displs);
    }

    free(recvbuf);

    MPI_Finalize();

    return 0;
}
