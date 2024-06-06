#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define INF 1e10
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

int main(int argc, char** argv) {
	clock_t start, end;
	double cpu_time_used;
	
	if(my_rank == 0){
		start = clock();
	}
	
	MPI_Init(NULL, NULL);
	int rank_size;
    MPI_Comm_size(MPI_COMM_WORLD, &rank_size);
	int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	if (argc < 2) {
        if (my_rank == 0) {
            printf("Incorrect usage. Expected usage: mpirun -np <p> ./floyd_warshall_mpi <inputMatrix sparse matrix>\n");
        }
        MPI_Finalize();
        return -1;
    }
	
	int V, E;
    FILE* ifile = fopen(argv[1], "r");
    if (ifile == NULL) {
        if (my_rank == 0) {
            printf("File not found.\n");
        }
        MPI_Finalize();
        return -1;
    }
	
	fscanf(ifile, "%d %d", &V, &E);
	
	// Allocate global data
    double** graph = (double**) malloc(V * sizeof(double*));
    for (int i = 0; i < V; i++) {
        graph[i] = (double*) malloc(V * sizeof(double));
        for (int j = 0; j < V; j++) {
            graph[i][j] = (i == j) ? 0 : INF;
        }
    }

    int source, target;
    double weight;
    // Store the weights in a matrix
    for (int i = 0; i < E; i++) {
        fscanf(ifile, "%d %d %lf", &source, &target, &weight);
        source--; target--; // Convert 1-based to 0-based index
        graph[source][target] = MIN(graph[source][target], weight);
    }
    fclose(ifile);
	
	int i, j, k, l, m;
	int count;
    for (k = 0; k < V; k++) {
        for (i = 0; i < V; i++) {
			count = 0;
			int indexes[V/rank_size];
			int changes[V/rank_size];
			for (j = (my_rank*(V/rank_size)); (j < (my_rank+1)*(V/rank_size)); j++) {
                if (graph[i][k] + graph[k][j] < graph[i][j]){
                    graph[i][j] = graph[i][k] + graph[k][j];
					indexes[count] = j;
					changes[count] = graph[i][j];
					count = count + 1;
				}
            }
			int counts[rank_size];
			MPI_Allgather(&count, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);
			
			int indexes_temp[V/rank_size];
			int changes_temp[V/rank_size];
			
			for(j = 0; j < rank_size; j++){
				if(counts[j] > 0){
					if(j == my_rank){
						for(l = 0; l < counts[j]; l++){
							indexes_temp[l] = indexes[l];
							changes_temp[l] = changes[l];
						}
					}
					MPI_Bcast(indexes_temp, counts[j], MPI_INT, j, MPI_COMM_WORLD);
					MPI_Bcast(changes_temp, counts[j], MPI_INT, j, MPI_COMM_WORLD);
					if(j != my_rank){
						for(l = 0; l < counts[j]; l++){
							graph[i][indexes_temp[l]] = changes_temp[l];
						}
					} 
				}
			}
        }
    }
	
	if(my_rank == 0){
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("time it took: %f \n", cpu_time_used);
	}
	
	if(my_rank == 0){
		FILE* outfile = fopen("outputMatrix.txt", "w");

		for (int i = 0; i < V; i++) {
			for (int j = 0; j < V; j++) {
				if (graph[i][j] == INF)
					fprintf(outfile,  "| %6s ", "INF");
				else
					fprintf(outfile, "| %6.2f ", graph[i][j]);
			}
			fprintf(outfile, "\n");
		}

		fclose(outfile);
	}
	
	MPI_Finalize();
    return 0;
}