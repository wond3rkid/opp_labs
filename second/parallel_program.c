#include "parallel_program.h"

void p_fill_initial_approximation(double *approximation, size_t N) {
    int i;
#pragma omp parallel for private(i)
    for (i = 0; i < N; i++) {
        approximation[i] = 0;
    }
}

void p_fill_matrix_vector(double **matrix, double *vector, size_t N) {
    int i, j;
#pragma omp parallel for shared(matrix, vector, N) private(i, j) default (none)
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            matrix[i][j] = (i == j) ? 2 : 1;
        }
        vector[i] = (double) N + 1;
    }
}

bool p_is_solved(const double **matrix, const double *vector, double *curr_approximation, size_t N) {
    double numerator_sqrt = p_get_vector_sqrt(vector, N);
    double *denominator = p_multiplication_matrix_vector(matrix, curr_approximation, N);
    p_subtracting_vectors(denominator, vector, N);
    double denominator_sqrt = p_get_vector_sqrt(denominator, N);
    free(denominator);
    return denominator_sqrt / numerator_sqrt < Epsilon;
}

double *p_multiplication_matrix_vector(const double **matrix, const double *vector, double *res, size_t N) {
    for (int i = 0; i < N; i++) {
        res[i] = 0;
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            res[i] += matrix[i][j] * vector[j];
        }
    }
}

void p_subtracting_vectors(double *curr, const double *vector, size_t N) {
    for (int i = 0; i < N; i++) {
        curr[i] -= vector[i];
    }
}

double *p_multiplication_tau_vector(const double *vector, size_t N) {

}

double p_get_vector_sqrt(const double *vector, size_t N);

double *p_get_next_x(const double **matrix, const double *vector, double *curr_approximation, size_t N);

void p_preparation_perfomance_free(size_t N);

void p_solve_equations(const double **matrix, const double *vector, double *initial_approximation, size_t N);

void p_print_result(double *res, size_t N);
