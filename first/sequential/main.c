#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

#define Epsilon 0.00001
#define Tau 0.00001
#define N 10000

double *create_matrix() {
    double *matrix = malloc(sizeof(double) * N * N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i * N + j] = (i == j) ? 2 : 1;
        }
    }
    return matrix;
}

double *create_vector() {
    double *vector = malloc(sizeof(double) * N);
    for (int i = 0; i < N; i++) {
        vector[i] = N + 1;
    }
    return vector;
}

void print_result(const double *vector) {
    printf("System of linear algebraic equations was solved. Check the result:\n");
    for (int i = 0; i < N; i++) {
        printf("%1.1f ", vector[i]);
    }
    printf("\n");
}

double get_vector_sqrt(const double *vector) {
    double res = 0;
    for (int i = 0; i < N; i++) {
        res += vector[i] * vector[i];
    }
    return pow(res, 0.5);
}

double *mult_tau_vector(const double *vector) {
    double *res = malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        res[i] = vector[i] * Tau;
    }
    return res;
}

double *mult_matrix_vector(const double *matrix, const double *vector) {
    double *result = calloc(N, sizeof(double));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            result[i] += matrix[i * N + j] * vector[j];
        }
    }
    return result;
}

double *subt_vectors(const double *curr, const double *vector) {
    double *res = calloc(N, sizeof(double));
    for (int i = 0; i < N; i++) {
        res[i] = curr[i] - vector[i];
    }
    return res;
}


bool is_solved(double *axb, double *b) {
    return get_vector_sqrt(axb) / get_vector_sqrt(b) < Epsilon;
}

void solve_equations(double *matrix, double *vector) {
    double *result = calloc(N, sizeof(double));
    while (true) {
        double *Ax = mult_matrix_vector(matrix, result);
        double *Axb = subt_vectors(Ax, vector);
        if (is_solved(Axb, vector)) {
            print_result(result);
            free(Ax);
            free(Axb);
            break;
        }
        double *tAxb = mult_tau_vector(Axb);
        result = subt_vectors(result, tAxb);
        free(tAxb);
    }
    free(result);
}

int main(int argc, char **argv) {
    double *vector = create_vector();
    double *matrix = create_matrix();
    clock_t start, end;
    start = clock();
    solve_equations(matrix, vector);
    end = clock();
    printf("Time taken : %ld sec\n", (end - start) / CLOCKS_PER_SEC);
    free(matrix);
    free(vector);
    return 0;
}