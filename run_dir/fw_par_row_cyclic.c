#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <float.h>
#include <time.h>

// Define infinity as a large number (used for initialization)
#define INF 1e10
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

void writeMatrixToFile(double** matrix, int V, const char* filename) {
    FILE* outfile = fopen(filename, "w");
    if (outfile == NULL) {
        fprintf(stderr, "Error opening output file!\n");
        return;
    }

    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            if (matrix[i][j] == INF)
                fprintf(outfile,  "| %6s ", "INF");
            else
                fprintf(outfile, "| %6.2f ", matrix[i][j]);
        }
        fprintf(outfile, "\n");
    }

    fclose(outfile);
}

int main(int argc, char** argv) {
    struct timespec start, end, floyd_start, floyd_end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int V, E;
    double** dist = NULL;

    if (rank == 0) {
        if (argc < 2) {
            printf("Incorrect usage. Expected usage: mpirun -np <p> ./floyd_warshall_mpi <inputMatrix sparse matrix>\n");
            MPI_Finalize();
            return -1;
        }

        FILE* ifile = fopen(argv[1], "r");
        if (ifile == NULL) {
            printf("File not found.\n");
            MPI_Finalize();
            return -1;
        }

        fscanf(ifile, "%d %d", &V, &E);

        if (V < size) {
            printf("Processes cannot be greater than vertices. Vertices: %d Comm size: %d\n", V, size);
            MPI_Finalize();
            return -1;
        }

        // Allocate global data
        dist = (double**) malloc(V * sizeof(double*));
        for (int i = 0; i < V; i++) {
            dist[i] = (double*) malloc(V * sizeof(double));
            for (int j = 0; j < V; j++) {
                dist[i][j] = (i == j) ? 0 : INF;
            }
        }

        int source, target;
        double weight;
        // Store the weights in a matrix
        for (int i = 0; i < E; i++) {
            fscanf(ifile, "%d %d %lf", &source, &target, &weight);
            source--; target--; // Convert 1-based to 0-based index
            dist[source][target] = MIN(dist[source][target], weight);
        }
        fclose(ifile);
    }

    // Broadcast the size of the matrix to all processes
    MPI_Bcast(&V, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int block_size = (V + size - 1) / size;

    // Allocate local data for each process
    double** local_dist = (double**) malloc(block_size * sizeof(double*));
    for (int i = 0; i < block_size; i++) {
        local_dist[i] = (double*) malloc(V * sizeof(double));
        for (int j = 0; j < V; j++) {
            local_dist[i][j] = INF;
        }
    }

    // Distribute rows in a block cyclic manner
    if (rank == 0) {
        for (int i = 0; i < V; i++) {
            int target_rank = i % size;
            if (target_rank == 0) {
                for (int j = 0; j < V; j++) {
                    local_dist[i / size][j] = dist[i][j];
                }
            } else {
                MPI_Send(dist[i], V, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD);
            }
        }
    } else {
        for (int i = 0; i < block_size; i++) {
            int global_row = i * size + rank;
            if (global_row < V) {
                MPI_Recv(local_dist[i], V, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &floyd_start);

    // Perform the local computation for this process's portion of the matrix
    for (int k = 0; k < V; k++) {
        // Calculate the source rank for the current k
        int root_rank = k % size;

        // Broadcast the k-th row to all processes
        double* k_row = (double*) malloc(V * sizeof(double));
        if ((k / size) < block_size && (k % size) == rank) {
            for (int j = 0; j < V; j++) {
                k_row[j] = local_dist[k / size][j];
            }
        }
        MPI_Bcast(k_row, V, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);

        for (int i = 0; i < block_size; i++) {
            int global_row = i * size + rank;
            if (global_row < V) {
                for (int j = 0; j < V; j++) {
                    local_dist[i][j] = MIN(local_dist[i][j], local_dist[i][k] + k_row[j]);
                }
            }
        }
        free(k_row);
    }

    // Gather the results back to the root process
    double* local_data = (double*) malloc(block_size * V * sizeof(double));
    for (int i = 0; i < block_size; i++) {
        for (int j = 0; j < V; j++) {
            local_data[i * V + j] = local_dist[i][j];
        }
    }

    double* global_data = NULL;
    if (rank == 0) {
        global_data = (double*) malloc(V * V * sizeof(double));
    }

    MPI_Gather(local_data, block_size * V, MPI_DOUBLE, global_data, block_size * V, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Convert the gathered 1D data into a 2D matrix
        clock_gettime(CLOCK_MONOTONIC, &floyd_end);

        double** global_dist = (double**) malloc(V * sizeof(double*));
        for (int i = 0; i < V; i++) {
            global_dist[i] = (double*) malloc(V * sizeof(double));
        }

        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                global_dist[i][j] = global_data[i * V + j];
            }
        }

        clock_gettime(CLOCK_MONOTONIC, &end);

        double total_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
        double floyd_time = (floyd_end.tv_sec - floyd_start.tv_sec) + (floyd_end.tv_nsec - floyd_start.tv_nsec) / 1e9;

        printf("Total Time measured: %.3f seconds.\n", total_time);
        printf("Floyd Time measured: %.3f seconds.\n", floyd_time);

        writeMatrixToFile(global_dist, V, "out_block_cyclic.txt");

        for (int i = 0; i < V; i++) {
            free(global_dist[i]);
        }
        free(global_dist);
    }

    MPI_Finalize();

    if (rank == 0) {
        for (int i = 0; i < V; i++) {
            free(dist[i]);
        }
        free(dist);
        free(global_data);
    }

    for (int i = 0; i < block_size; i++) {
        free(local_dist[i]);
    }
    free(local_dist);
    free(local_data);

    return 0;
}
