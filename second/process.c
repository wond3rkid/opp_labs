#include "math.h"
#include "process.h"

void fill_vector_initial_approximation(double *approximation, size_t N) {
    for (int i = 0; i < N; i++) {
        approximation[i] = 10;
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

void subtracting_vectors(double *curr, const double *vector, size_t N) {
    for (int i = 0; i < N; i++) {
        curr[i] -= vector[i];
    }
}

double *multiplication_tau_vector(const const double *vector, size_t N) {
    double *result = malloc(sizeof(result) * N);
    for (int i = 0; i < N; i++) {
        result[i] = vector[i] * Tau;
    }
    return result;
}

double get_vector_sqrt(const double *vector, size_t N) {
    double ans = 0;
    for (int i = 0; i < N; i++) {
        ans += vector[i] * vector[i];
    }
    return (double) pow(ans, 0.5);
}

bool is_solved(const double **matrix, const double *vector, double *curr_approximation, size_t N) {
    double numerator_sqrt = get_vector_sqrt(vector, N);
    double *denominator = multiplication_matrix_vector(matrix, curr_approximation, N);
    subtracting_vectors(denominator, vector, N);
    print_result(denominator, N);
    double denominator_sqrt = get_vector_sqrt(denominator, N);
    print_result(vector, N);
    free(denominator);
    return denominator_sqrt / numerator_sqrt < Epsilon;
}

double *get_next_x(const double **matrix, const double *vector, double *curr_approximation, size_t N) {
    double *tmp_vect = multiplication_matrix_vector(matrix, curr_approximation, N);
    subtracting_vectors(tmp_vect, vector, N);
    double *tmp_curr = multiplication_tau_vector(tmp_vect, N);
    subtracting_vectors(curr_approximation, tmp_curr, N);
    return curr_approximation;
}

void preparation_perfomance_free(size_t N) {
    double **matrix = malloc(sizeof(*matrix) * N);
    for (int i = 0; i < N; i++) {
        matrix[i] = malloc(sizeof(matrix[i]) * N);
    }
    double *vector = malloc(sizeof(vector) * N);
    fill_matrix_vector(matrix, vector, N);
    double *initial_approximation = malloc(sizeof(vector) * N);
    fill_vector_initial_approximation(initial_approximation, N);
    solve_equations(matrix, vector, initial_approximation, N);
    print_result(initial_approximation, N);
    for (int i = 0; i < N; i++) {
        free(matrix[i]);
    }
    free(matrix);
    free(vector);
}

void solve_equations(const double **matrix, const double *vector, double *initial_approximation, size_t N) {
    do {
        get_next_x(matrix, vector, initial_approximation, N);
        printf("\n");
    } while (!is_solved(matrix, vector, initial_approximation, N));
}

void print_result(double *result, size_t N) {
    for (size_t i = 0; i < N; i++) {
        printf("res[%d] = %f\n", i, result[i]);
    }
}
