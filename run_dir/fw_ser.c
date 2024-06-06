#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdbool.h>
#include <time.h>

#define INF DBL_MAX

typedef struct {
    double weight;
    bool infinite;
} Distance;

// Function prototypes
void setWeight(Distance* d, double weight);
char* getWeight(Distance* d, char* buffer);
bool isInfinite(Distance* d);
bool isZero(Distance* d);
bool compareDistance(Distance* d1, Distance* d2);
Distance addDistance(Distance* d1, Distance* d2);
char* obtainPath(int i, int j, int **parent, Distance **dist, char* buffer);

void setWeight(Distance* d, double weight) {
    d->weight = weight;
    d->infinite = false;
}

char* getWeight(Distance* d, char* buffer) {
    if (d->infinite) 
        sprintf(buffer, "inf");
    else
        sprintf(buffer, "%.2f", d->weight);
    return buffer;
}

bool isInfinite(Distance* d) {
    return d->infinite;
}

bool isZero(Distance* d) {
    return (!d->infinite && d->weight == 0);
}

bool compareDistance(Distance* d1, Distance* d2) {
    if (d2->infinite)
        return false;
    else if (d1->infinite)
        return true;
    else if (d1->weight > d2->weight)
        return true;
    return false;
}

Distance addDistance(Distance* d1, Distance* d2) {
    Distance d;
    d.infinite = true;
    if (d1->infinite || d2->infinite)
        return d;
    d.weight = d1->weight + d2->weight;
    d.infinite = false;
    return d;
}

char* obtainPath(int i, int j, int **parent, Distance **dist, char* buffer) {
    static char result[1024];
    if (dist[i][j].infinite)
        sprintf(result, " no path to ");
    else if (parent[i][j] == i)
        sprintf(result, " ");
    else {
        sprintf(result, "%s%d%s", obtainPath(i, parent[i][j], parent, dist, buffer), parent[i][j] + 1, obtainPath(parent[i][j], j, parent, dist, buffer));
    }
    return result;
}

int main(int argc, char** argv) {
    clock_t start = clock();
    if (argc < 2) {
        printf("Check README for usage.\n");
        exit(-1);
    }

    FILE* ifile = fopen(argv[1], "r");
    if (!ifile) {
        printf("File not found.\n");
        exit(-1);
    }

    // number of vertices and edges
    int V, E;
    fscanf(ifile, "%d %d", &V, &E);

    // Matrices declared and initialised to infinity and zero respectively
    Distance **dist = (Distance**)malloc(V * sizeof(Distance*));
    for (int i = 0; i < V; i++) {
        dist[i] = (Distance*)malloc(V * sizeof(Distance));
        for (int j = 0; j < V; j++) {
            dist[i][j].weight = 0;
            dist[i][j].infinite = true;
        }
    }

    int **parent = (int**)malloc(V * sizeof(int*));
    for (int i = 0; i < V; i++) {
        parent[i] = (int*)malloc(V * sizeof(int));
    }

    // Read edges from input file
    for (int i = 0; i < E; i++) {
        int u, v;
        double w;
        fscanf(ifile, "%d %d %lf", &u, &v, &w);
        setWeight(&dist[u - 1][v - 1], w);
        parent[u - 1][v - 1] = u - 1;
    }
    fclose(ifile);

    // Path from vertex to itself is set
    for (int i = 0; i < V; i++) {
        setWeight(&dist[i][i], 0);
        parent[i][i] = i;
    }

    clock_t floyd_start = clock();
    // Actual Floyd Warshall Algorithm
    for (int k = 0; k < V; k++) {
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                Distance temp = addDistance(&dist[i][k], &dist[k][j]);
                if (compareDistance(&dist[i][j], &temp)) {
                    dist[i][j] = temp;
                    parent[i][j] = parent[k][j];
                }
            }
        }
    }

    // Check for negative cycles
    for (int i = 0; i < V; i++) {
        if (!isZero(&dist[i][i])) {
            printf("Negative cycle at : %d\n", i + 1);
            return 0;
        }
    }

    clock_t floyd_end = clock();

    FILE* outfile = fopen("out_ser.txt", "w");
    if (!outfile) {
        fprintf(stderr, "Error opening output file!\n");
        return 1; // Exit with an error code
    }

    // Print all paths in the specified format
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            char buffer[20];
            if (j > 0) {
                fprintf(outfile, " ");
            }
            fprintf(outfile, "| %6s", getWeight(&dist[i][j], buffer));
        }
        fprintf(outfile, "|\n");
    }
    fclose(outfile);

    for (int i = 0; i < V; i++) {
        free(dist[i]);
        free(parent[i]);
    }
    free(dist);
    free(parent);

    clock_t end = clock();
    double total_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    double floyd_time = ((double)(floyd_end - floyd_start)) / CLOCKS_PER_SEC;

    printf("Total Time measured: %.3f seconds.\n", total_time);
    printf("Floyd Time measured: %.3f seconds.\n", floyd_time);

    return 0;
}
