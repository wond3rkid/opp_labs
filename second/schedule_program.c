#include "schedule_program.h"
#include "parallel_for_program.h"

void sc_fill_initial_approximation(double *approximation, size_t N) {
#pragma omp parallel for schedule (static, 500)
    for (int i = 0; i < N; i++) {
        approximation[i] = 0;
    }
}

void sc_fill_matrix_vector(double **matrix, double *vector, size_t N) {
#pragma omp parallel for schedule (static, 500)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j] = (i == j) ? 2 : 1;
        }
        vector[i] = (double) N + 1;
    }
}

bool sc_is_solved(const double **matrix, const double *vector, double *curr_approximation, size_t N) {
    double numerator_sqrt = sc_get_vector_sqrt(vector, N);
    double *denominator = malloc(sizeof(denominator) * N);
    sc_multiplication_matrix_vector(matrix, curr_approximation, denominator, N);
    sc_subtracting_vectors(denominator, vector, N);
    double denominator_sqrt = sc_get_vector_sqrt(denominator, N);
    free(denominator);
    return denominator_sqrt / numerator_sqrt < Epsilon;
}

void sc_multiplication_matrix_vector(const double **matrix, const double *vector, double *res, size_t N) {
#pragma omp parallel for schedule (static, 500)
    for (int i = 0; i < N; i++) {
        res[i] = 0;
    }
#pragma omp parallel for schedule (static, 500)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            res[i] += matrix[i][j] * vector[j];
        }
    }
}

void sc_subtracting_vectors(double *curr, const double *vector, size_t N) {
#pragma omp parallel for schedule (static, 500)
    for (int i = 0; i < N; i++) {
        curr[i] -= vector[i];
    }
}

void sc_multiplication_tau_vector(const double *vector, double *result, size_t N) {
#pragma omp parallel for schedule (static, 500)
    for (int i = 0; i < N; i++) {
        result[i] = vector[i] * Tau;
    }
}

double sc_get_vector_sqrt(const double *vector, size_t N) {
    double ans = 0;
#pragma omp parallel for reduction(+: ans)
    for (int i = 0; i < N; i++) {
        ans += vector[i] * vector[i];
    }
    return (double) pow(ans, 0.5);
}


void sc_preparation_perfomance_free(size_t N) {
    double **matrix = malloc(sizeof(*matrix) * N);
    for (int i = 0; i < N; i++) {
        matrix[i] = malloc(sizeof(matrix[i]) * N);
    }
    double *vector = malloc(sizeof(vector) * N);
    sc_fill_matrix_vector(matrix, vector, N);
    double *initial_approximation = malloc(sizeof(vector) * N);
    sc_fill_initial_approximation(initial_approximation, N);
    sc_solve_equations((const double **) matrix, vector, initial_approximation, N);
   // sc_print_result(initial_approximation, N);
    for (int i = 0; i < N; i++) {
        free(matrix[i]);
    }
    free(initial_approximation);
    free(matrix);
    free(vector);
}

void sc_solve_equations(const double **matrix, const double *vector, double *initial_approximation, size_t N) {
    do {
        double *tmp_vect = malloc(sizeof(tmp_vect) * N);
        sc_multiplication_matrix_vector(matrix, initial_approximation, tmp_vect, N);
        sc_subtracting_vectors(tmp_vect, vector, N);
        double *tmp_curr = malloc(sizeof(tmp_curr) * N);
        sc_multiplication_tau_vector(tmp_vect, tmp_curr, N);
        sc_subtracting_vectors(initial_approximation, tmp_curr, N);
        free(tmp_curr);
        free(tmp_vect);
    } while (!sc_is_solved(matrix, vector, initial_approximation, N));

}


void sc_print_result(double *res, size_t N) {
    for (int i = 0; i < N; i++) {
        printf("res[%d] = %f\n", i, res[i]);
    }
}