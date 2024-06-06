#include <stdio.h>
#include <stdlib.h>
#include <float.h>

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

    // Perform Floyd-Warshall algorithm
    floydWarshall(n, dist);

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

    return 0;
}
