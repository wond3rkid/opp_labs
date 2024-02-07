#include <stdio.h>


void fillMatrix(double **matrix_A, double *matrix_B, size_t N) {
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            matrix_A[i][j] = (i == j) ? 2 : 1;
        } matrix_B[i] = N + 1;
    }
}

void iterationMethod();

int main(int argc, char **argv) {
    printf("Hello, World!\n");
    return 0;
}
