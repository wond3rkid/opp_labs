#include <stdio.h>
#include <malloc.h>
#include <assert.h>

#define N 4 //n1
#define M 5 //n2
#define P 3 //n3

double *init_matrix(int r, int c) {
    double *matrix = calloc(r * c, sizeof(double));
    assert(matrix != NULL);
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            matrix[i * c + j] = (j == i ? 2 : 1);
        }
    }
    return matrix;
}


void print_matrix(int r, int c, const double *matrix) {
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            printf("%1.2f \t", matrix[i * c + j]);
        }
        printf("\n");
    }
}

double *multiplicate_matrix(const double *A, const double *B, int n1, int n2, int n3) {
    double *result = calloc(n1 * n3, sizeof(double));
    assert(result != NULL);
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n3; j++) {
            result[i * n3 + j] = 0;
            for (int k = 0; k < n2; k++) {
                result[i * n1 + j] += A[i * n1 + k] * B[k * n3 + j];
            }
        }
    }
    return result;
}

void do_task() {
    double *matrix_a = init_matrix(N, M);
    double *matrix_b = init_matrix(M, P);
    print_matrix(N, M, matrix_a);
    printf("\n");
    print_matrix(M, P, matrix_b);
    double *res = multiplicate_matrix(matrix_a, matrix_b, N, M, P);
    printf("\n");

    print_matrix(N, P, res);
    free(res);
    free(matrix_b);
    free(matrix_a);
}

int main() {
    do_task();
    return 0;
}