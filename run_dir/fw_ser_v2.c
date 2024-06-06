#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#define MAX_NODES 4000

void floydWarshall(int n, double dist[][MAX_NODES]) {
    for (int k = 1; k <= n; k++) {
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
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
    if (fscanf(file, "%d %d", &n, &m) != 2) {
        fprintf(stderr, "Error: failed to read number of nodes and edges.\n");
        fclose(file);
        return 1;
    }

    if (n > MAX_NODES || n <= 0) {
        fprintf(stderr, "Error: number of nodes is out of bounds. n=%d\n", n);
        fclose(file);
        return 1;
    }

    double (*dist)[MAX_NODES] = (double (*)[MAX_NODES]) malloc(sizeof(double[MAX_NODES][MAX_NODES]));
    if (dist == NULL) {
        fprintf(stderr, "Error: failed to allocate memory for distance matrix.\n");
        fclose(file);
        return 1;
    }

    // Initialize distance matrix
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
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
        if (fscanf(file, "%d %d %lf", &u, &v, &weight) != 3) {
            fprintf(stderr, "Error: failed to read edge %d.\n", i + 1);
            free(dist);
            fclose(file);
            return 1;
        }
        if (u < 1 || u > n || v < 1 || v > n) {
            fprintf(stderr, "Error: edge %d has out of bounds vertices. u=%d, v=%d, n=%d\n", i + 1, u, v, n);
            free(dist);
            fclose(file);
            return 1;
        }
        dist[u][v] = weight;
    }

    fclose(file);

    // Perform Floyd-Warshall algorithm
    floydWarshall(n, dist);

    // Output the result to out_ser.txt
    FILE *outfile = fopen("out_ser.txt", "w");
    if (outfile == NULL) {
        fprintf(stderr, "Error: could not open output file.\n");
        free(dist);
        return 1;
    }

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (dist[i][j] == DBL_MAX) {
                fprintf(outfile, "INF ");
            } else {
                fprintf(outfile, "%.6lf ", dist[i][j]);
            }
        }
        fprintf(outfile, "\n");
    }

    fclose(outfile);
    free(dist);

    return 0;
}
