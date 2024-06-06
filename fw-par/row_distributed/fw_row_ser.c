#include <stdio.h>
#include <stdlib.h>

void floyd_warshall_serial(int num_of_vertices, double *matrix) {
    int i, j, k;
    for (k = 0; k < num_of_vertices; k++) {
        for (i = 0; i < num_of_vertices; i++) {
            for (j = 0; j < num_of_vertices; j++) {
                if (matrix[i * num_of_vertices + j] > matrix[i * num_of_vertices + k] + matrix[k * num_of_vertices + j]) {
                    matrix[i * num_of_vertices + j] = matrix[i * num_of_vertices + k] + matrix[k * num_of_vertices + j];
                }
            }
        }
    }
}

int main() {
    int num_of_vertices = 12;
    double matrix[12 * 12] = {
        0, 4, 2, 1, 7, 6, 5, 1, 7, 2, 9, 3,
        5, 0, 3, 2, 0, 8, 1, 9, 1, 8, 1, 8,
        3, 0, 0, 6, 2, 8, 1, 6, 3, 4, 7, 0,
        0, 2, 2, 0, 4, 2, 1, 0, 4, 4, 3, 5,
        3, 4, 4, 5, 0, 6, 3, 6, 6, 2, 2, 9,
        1, 4, 5, 4, 8, 0, 5, 9, 4, 7, 7, 9,
        9, 9, 0, 4, 4, 3, 0, 7, 8, 4, 3, 1,
        1, 6, 7, 7, 9, 9, 6, 0, 3, 1, 5, 2,
        4, 1, 1, 9, 9, 9, 9, 9, 0, 9, 4, 3,
        3, 4, 0, 1, 9, 3, 2, 0, 0, 0, 7, 9,
        9, 4, 0, 3, 5, 6, 5, 0, 7, 7, 0, 7,
        6, 8, 6, 5, 8, 1, 8, 1, 5, 9, 2, 0
    };

    floyd_warshall_serial(num_of_vertices, matrix);

    for (int i = 0; i < num_of_vertices; i++) {
        for (int j = 0; j < num_of_vertices; j++) {
            printf("| %lf ", matrix[i * num_of_vertices + j]);
        }
        printf("|\n");
    }

    return 0;
}
