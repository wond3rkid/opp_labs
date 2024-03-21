#include "parallel_section_program.h"

void s_fill_initial_approximation(double *approximation, size_t N) {
    for (int i = 0; i < N; i++) {
        approximation[i] = 0;
    }
}

void s_fill_matrix_vector(double **matrix, double *vector, size_t N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j] = (i == j) ? 2 : 1;
        }
        vector[i] = (double) N + 1;
    }
}

bool s_is_solved(const double **matrix, const double *vector, double *curr_approximation, size_t N) {
    double numerator_sqrt = s_get_vector_sqrt(vector, N);
    double *denominator = malloc(sizeof(denominator) * N);
    //fprintf(stderr, "HERE");
    s_multiplication_matrix_vector(matrix, curr_approximation, denominator, N);
    s_subtracting_vectors(denominator, vector, N);
    double denominator_sqrt = s_get_vector_sqrt(denominator, N);
    free(denominator);
    // there is 1 thread
    // error is here

    //fprintf(stderr, "%f %f \n", denominator_sqrt, denominator_sqrt / numerator_sqrt );
    return denominator_sqrt / numerator_sqrt < Epsilon;
}

void s_multiplication_matrix_vector(const double **matrix, const double *vector, double *res, size_t N) {
// there is 12 threads
#pragma omp for
    for (int i = 0; i < N; i++) {
/*        if (i < 15) {
            fprintf(stderr, "s %d - %d thread\n", i, omp_get_thread_num());
        }*/
        for (int j = 0; j < N; j++) {
            res[i] += matrix[i][j] * vector[j];
        }
    }
}

void s_subtracting_vectors(double *curr, const double *vector, size_t N) {
    // there is 12 threads
#pragma omp for
    for (int i = 0; i < N; i++) {
        curr[i] -= vector[i];
    }
}

void s_multiplication_tau_vector(const double *vector, double *result, size_t N) {
    // there is 12 threads
#pragma omp for
    for (int i = 0; i < N; i++) {
        result[i] = vector[i] * Tau;
    }
}

double s_get_vector_sqrt(const double *vector, size_t N) {
    double ans = 0;
    for (int i = 0; i < N; i++) {
        ans += vector[i] * vector[i];
    }
    return (double) pow(ans, 0.5);
}

void s_preparation_perfomance_free(size_t N) {
    double **matrix = malloc(sizeof(*matrix) * N);
    for (int i = 0; i < N; i++) {
        matrix[i] = malloc(sizeof(matrix[i]) * N);
    }
    double *vector = malloc(sizeof(vector) * N);

    s_fill_matrix_vector(matrix, vector, N);

    double *initial_approximation = malloc(sizeof(vector) * N);

    s_fill_initial_approximation(initial_approximation, N);

    s_solve_equations((const double **) matrix, vector, initial_approximation, N);

    //s_print_result(initial_approximation, N);

    for (int i = 0; i < N; i++) {
        free(matrix[i]);
    }
    free(initial_approximation);
    free(matrix);
    free(vector);
}

void s_solve_equations(const double **matrix, const double *vector, double *initial_approximation, size_t N) {
    do {
        double *tmp_vect = malloc(sizeof(tmp_vect) * N);
        double *tmp_curr = malloc(sizeof(tmp_curr) * N);
        for (int i = 0; i < N; i++) {
            tmp_vect[i] = 0;
        }
#pragma omp parallel
        {
            s_multiplication_matrix_vector(matrix, initial_approximation, tmp_vect, N);
#pragma omp barrier
            s_subtracting_vectors(tmp_vect, vector, N);
#pragma omp barrier
            s_multiplication_tau_vector(tmp_vect, tmp_curr, N);
#pragma omp barrier
            s_subtracting_vectors(initial_approximation, tmp_curr, N);
        }
        free(tmp_curr);
        free(tmp_vect);
    } while (!s_is_solved(matrix, vector, initial_approximation, N));
}

void s_print_result(double *res, size_t N) {
    for (int i = 0; i < N; i++) {
        printf("res[%d] = %f\n", i + 1, res[i]);
    }
}
