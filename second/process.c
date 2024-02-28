
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

double *multiplication_matrix_vector(const double **matrix, double *curr_approximation, size_t N) {
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
        curr[i] -= vector[i];
    }
}

double *multiplication_tau_vector(const double *vector, size_t N) {
    double *result = malloc(sizeof(result) * N);
    for (int i = 0; i < N; i++) {
        result[i] = vector[i] * Tau;
    }
    return result;
}

double get_vector_sqrt(double *vector, size_t N) {
    double ans = 0;
    for (int i = 0; i < N; i++) {
        ans += vector[i] * vector[i];
    }
    return powl(ans, 0.5);
}

bool is_solved(const double **matrix, double *vector, double *curr_approximation, size_t N) {
    double numerator_sqrt = get_vector_sqrt(vector, N);
    double *denominator = multiplication_matrix_vector(matrix, curr_approximation, N);
    subtracting_vectors(denominator, vector, N);
    double denominator_sqrt = get_vector_sqrt(denominator, N);
    free(denominator);
    return numerator_sqrt / denominator_sqrt < Epsilon;
}

double *iterationMethod(const double **matrix_A, double *matrix_B, size_t N) {

}

double *get_next_x(const double **matrix, const double *vector, double *curr_approximation, size_t N) {
    double *tmp_vect = multiplication_matrix_vector(matrix, curr_approximation, N);
    subtracting_vectors(tmp_vect, vector, N);
    double *tmp_curr = multiplication_tau_vector(tmp_vect, N);
    subtracting_vectors(curr_approximation, tmp_curr, N);
    return curr_approximation;
}

void do_algorithm(size_t N) {



}


void print_result(double *result, size_t N) {
    for (size_t i = 0; i < N; i++) {
        printf("res[i] = %d\n", result[i]);
    }
}