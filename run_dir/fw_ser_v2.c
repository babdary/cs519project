#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>

void floydWarshall(int n, double **dist) {
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (dist[i][k] < DBL_MAX && dist[k][j] < DBL_MAX && dist[i][k] + dist[k][j] < dist[i][j]) {
                    dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
    }
}

int main() {
    FILE *file = fopen("inputMatrix.txt", "r");
    if (file == NULL) {
        fprintf(stderr, "Error: could not open input file.\n");
        return 1;
    }

    int n, m;
    fscanf(file, "%d %d", &n, &m);

    // Allocate memory for the distance matrix
    double **dist = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        dist[i] = (double *)malloc(n * sizeof(double));
    }

    // Initialize distance matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                dist[i][j] = 0.0;
            } else {
                dist[i][j] = DBL_MAX;
            }
        }
    }

    // Read the edges
    for (int i = 0; i < m; i++) {
        int u, v;
        double weight;
        fscanf(file, "%d %d %lf", &u, &v, &weight);
        dist[u - 1][v - 1] = weight; // Adjust indices to start from 0
    }

    fclose(file);

    // Measure the start time
    clock_t start_time = clock();

    // Perform Floyd-Warshall algorithm
    floydWarshall(n, dist);

    // Measure the end time
    clock_t end_time = clock();

    // Output the result to out_ser.txt
    FILE *outfile = fopen("out_ser.txt", "w");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (dist[i][j] == DBL_MAX) {
                fprintf(outfile, "INF ");
            } else {
                fprintf(outfile, "%.6lf ", dist[i][j]);
            }
        }
        fprintf(outfile, "\n");
    }
    fclose(outfile);

    // Free allocated memory
    for (int i = 0; i < n; i++) {
        free(dist[i]);
    }
    free(dist);

    // Calculate and print the total runtime
    double total_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Total Time measured: %.3lf seconds.\n", total_time);

    return 0;
}
