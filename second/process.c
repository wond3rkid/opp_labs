
#include "process.h"

void fill_vector_initial_approximation(double *approximation, size_t N) {
    for (int i = 0; i < N; i++) {
        approximation[i] = 1;
    }
}

void fill_matrix_vector(double **matrix, double *vector, size_t N) {
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            matrix[i][j] = (i == j) ? 2 : 1;
        }
        vector[i] = (double) N + 1;
    }
}

double *multiplication_matrix_vector(double **matrix, double *curr_approximation, size_t N) {
    double *res = malloc(sizeof(res) * N);
    for (int i = 0; i < N; i++) {
        res[i] = 0;
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            res[i] += matrix[i][j] * curr_approximation[j];
        }
    }
    return res;
}

void subtracting_vectors(double *curr, double *vector, size_t N) {
    for (int i = 0; i < N; i++) {
        curr[i] += vector[i];
    }
}

double **multiplication_tau_matrix(double **matrix, size_t N) {
    double **result = malloc(sizeof(*result) * N);
    for (int i = 0; i < N; i++) {
        result[i] = malloc(sizeof(result) * N);
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            result[i][j] = matrix[i][j] * Tau;
        }
    }
}

double get_vector_sqrt(double *vector, size_t N) {
    double ans = 0;
    for (int i = 0; i < N; i++) {
        ans += vector[i] * vector[i];
    }
    return powl(ans, 0.5);
}

bool is_solved(double **matrix, double *vector, double *curr_approximation, size_t N) {
    double numerator_sqrt = get_vector_sqrt(vector, N);
    double *denominator = multiplication_matrix_vector(matrix, curr_approximation, N);
    subtracting_vectors(denominator, vector, N);
    double denominator_sqrt = get_vector_sqrt(denominator, N);
    free(denominator);
    return numerator_sqrt / denominator_sqrt < Epsilon;
}

double *iterationMethod(double **matrix_A, double *matrix_B, size_t N) {

}

void do_algorithm(size_t N) {
    double **matrixA = malloc(sizeof(*matrixA) * N);
    double *matrixB = malloc(sizeof(matrixB) * N);

}


void print_result(double *res, size_t N) {
    for (size_t i = 0; i < N; i++) {
        printf("res[i] = %d\n", res[i]);
    }
}