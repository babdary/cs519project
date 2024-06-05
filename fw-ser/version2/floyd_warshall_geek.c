#include <stdio.h>
#include <stdlib.h>

#define INF 99999.0

void printSolution(double **dist, int V);

void floydWarshall(double **dist, int V) {
    int i, j, k;
    for (k = 0; k < V; k++) {
        for (i = 0; i < V; i++) {
            for (j = 0; j < V; j++) {
                if (dist[i][k] + dist[k][j] < dist[i][j])
                    dist[i][j] = dist[i][k] + dist[k][j];
            }
        }
    }
    printSolution(dist, V);
}

void printSolution(double **dist, int V) {
    printf("The following matrix shows the shortest distances between every pair of vertices \n");
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            if (dist[i][j] == INF)
                printf("%7s ", "INF");
            else
                printf("%d ", (int)dist[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Usage: %s <input file>\n", argv[0]);
        return 1;
    }

    FILE *file = fopen(argv[1], "r");
    if (file == NULL) {
        printf("Error opening file\n");
        return 1;
    }

    int V;
    fscanf(file, "%d", &V);

    double **graph = (double **)malloc(V * sizeof(double *));
    for (int i = 0; i < V; i++) {
        graph[i] = (double *)malloc(V * sizeof(double));
        for (int j = 0; j < V; j++) {
            if (i == j)
                graph[i][j] = 0.0;
            else
                graph[i][j] = INF;
        }
    }

    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            fscanf(file, "%lf", &graph[i][j]);
            if (i != j && graph[i][j] == 0) {
                graph[i][j] = INF; // Handle the case where no edge exists
            }
        }
    }

    fclose(file);

    floydWarshall(graph, V);

    for (int i = 0; i < V; i++) {
        free(graph[i]);
    }
    free(graph);

    return 0;
}
