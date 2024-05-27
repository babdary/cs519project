#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

typedef struct {
    int source;
    int target;
    int weight;
} Edge;

void generate_graph(Edge *graph, int nodes, int edges) {
    int i = 0;
    while (i < edges) {
        int source = rand() % nodes + 1;
        int target = rand() % nodes + 1;
        int weight = rand() % 100;

        if (source != target) {
            bool duplicate = false;
            for (int j = 0; j < i; j++) {
                if (graph[j].source == source && graph[j].target == target) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                graph[i].source = source;
                graph[i].target = target;
                graph[i].weight = weight;
                i++;
            }
        }
    }
}

void save_graph(const char *filename, Edge *graph, int nodes, int edges) {
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        perror("Failed to open file");
        exit(EXIT_FAILURE);
    }

    fprintf(f, "%d\n%d\n", nodes, edges);
    for (int i = 0; i < edges; i++) {
        fprintf(f, "%d %d %d\n", graph[i].source, graph[i].target, graph[i].weight);
    }

    fclose(f);
}

void save_adjacency_matrix(const char *filename, Edge *graph, int nodes, int edges) {
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        perror("Failed to open file");
        exit(EXIT_FAILURE);
    }

    int **matrix = (int **)malloc(nodes * sizeof(int *));
    for (int i = 0; i < nodes; i++) {
        matrix[i] = (int *)calloc(nodes, sizeof(int));
    }

    for (int i = 0; i < edges; i++) {
        int source = graph[i].source - 1;
        int target = graph[i].target - 1;
        int weight = graph[i].weight;
        matrix[source][target] = weight;
    }

    fprintf(f, "%d\n", nodes);
    for (int i = 0; i < nodes; i++) {
        for (int j = 0; j < nodes; j++) {
            fprintf(f, "%d ", matrix[i][j]);
        }
        fprintf(f, "\n");
    }

    for (int i = 0; i < nodes; i++) {
        free(matrix[i]);
    }
    free(matrix);

    fclose(f);
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <nodes> <edges> [generate_matrix]\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    int nodes = atoi(argv[1]);
    int edges = atoi(argv[2]);
    bool generate_matrix = argc > 3 ? atoi(argv[3]) : false;

    srand(time(NULL));

    Edge *graph = (Edge *)malloc(edges * sizeof(Edge));
    if (graph == NULL) {
        perror("Failed to allocate memory");
        exit(EXIT_FAILURE);
    }

    generate_graph(graph, nodes, edges);
    save_graph("inputMatrix.txt", graph, nodes, edges);

    if (generate_matrix) {
        save_adjacency_matrix("adjacencyMatrix.txt", graph, nodes, edges);
    }

    free(graph);
    return 0;
}
