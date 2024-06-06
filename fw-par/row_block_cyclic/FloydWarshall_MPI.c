#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <float.h>
#include <time.h>

// Define infinity as a large number (used for initialization)
#define INF 1e10
//#define min(a,b) if(a<b) ? a:b
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

    if (argc < 2) {
        if (rank == 0) {
            printf("Incorrect usage. Expected usage: mpirun -np <p> ./floyd_warshall_mpi <inputMatrix sparse matrix>\n");
        }
        MPI_Finalize();
        return -1;
    }

    int V, E;
    FILE* ifile = fopen(argv[1], "r");
    if (ifile == NULL) {
        if (rank == 0) {
            printf("File not found.\n");
        }
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
    double** dist = (double**) malloc(V * sizeof(double*));
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

    int block_size = (V + size - 1) / size;  // Calculate block size +size -1 is used to handle boundary conditions

    clock_gettime(CLOCK_MONOTONIC, &floyd_start);

    // Allocate local data for each process
    double** local_dist = (double**) malloc(block_size * sizeof(double*));
    for (int i = 0; i < block_size; i++) {
        local_dist[i] = (double*) malloc(V * sizeof(double));
        for (int j = 0; j < V; j++) {
            local_dist[i][j] = INF;
        }
    }

    double* local_data = (double*) malloc(block_size * (V + 1) * sizeof(double));
    int cc = 0;
    // Distribute rows in a block cyclic manner
    for (int i = 0; i < block_size; i++) {
        int global_row = i * size + rank;
        if (global_row < V) {
            for (int j = 0; j < V; j++) {
                local_dist[i][j] = dist[global_row][j];
                local_data[i * V + j + cc] = dist[global_row][j];
            }
            local_data[V + i * V + cc] = global_row;
            cc++;
        }
    }

    // Perform the local computation for this process's portion of the matrix
    for (int k = 0; k < V; k++) {
        // Calculate the source rank for the current k
        int root_rank = k % size;

        // Broadcast the k-th row to all processes
        double* k_row = (double*) malloc(V * sizeof(double));
        for (int j = 0; j < V; j++) {
            k_row[j] = INF;
        }
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
                    if (local_dist[i][j] > local_dist[i][k] + k_row[j]) {
                        local_dist[i][j] = local_dist[i][k] + k_row[j];
                    }
                }
            }
        }
        free(k_row);
    }

    // Update local_data from local_dist before gathering
    cc = 0;
    for (int i = 0; i < block_size; i++) {
        int global_row = i * size + rank;
        if (global_row < V) {
            for (int j = 0; j < V; j++) {
                local_data[i * V + j + cc] = local_dist[i][j];
            }
            cc++;
        }
    }

    // Prepare the buffer for gathering the results
    double* global_data = NULL;
    if (rank == 0) {
        global_data = (double*) malloc(V * (V + 1) * sizeof(double));
        for (int i = 0; i < V * (V + 1); i++) {
            global_data[i] = INF;
        }
    }

    // Gather the results back to the root process
    MPI_Gather(local_data, block_size * (V + 1), MPI_DOUBLE,
               global_data, block_size * (V + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Convert the gathered 1D data into a 2D matrix
        cc = 0;
        clock_gettime(CLOCK_MONOTONIC, &floyd_end);

        double** global_dist = (double**) malloc(V * sizeof(double*));
        for (int i = 0; i < V; i++) {
            global_dist[i] = (double*) malloc(V * sizeof(double));
            for (int j = 0; j < V; j++) {
                global_dist[i][j] = INF;
            }
        }

        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                global_dist[(int)global_data[i * V + V + cc]][j] = global_data[i * V + j + cc];
            }
            cc++;
        }

        clock_gettime(CLOCK_MONOTONIC, &end);

        double total_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
        double floyd_time = (floyd_end.tv_sec - floyd_start.tv_sec) + (floyd_end.tv_nsec - floyd_start.tv_nsec) / 1e9;

        printf("Total Time measured: %.3f seconds.\n", total_time);
        printf("Floyd Time measured: %.3f seconds.\n", floyd_time);

        writeMatrixToFile(global_dist, V, "output_matrix.txt");

        for (int i = 0; i < V; i++) {
            free(global_dist[i]);
        }
        free(global_dist);
    }

    MPI_Finalize();

    for (int i = 0; i < V; i++) {
        free(dist[i]);
    }
    free(dist);

    for (int i = 0; i < block_size; i++) {
        free(local_dist[i]);
    }
    free(local_dist);
    free(local_data);
    if (rank == 0) {
        free(global_data);
    }

    return 0;
}
