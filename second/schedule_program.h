#ifndef SECOND_SCHEDULE_PROGRAM_H
#define SECOND_SCHEDULE_PROGRAM_H

#include <stdio.h>
#include <stdbool.h>

#define Epsilon 0.00001
#define Tau 0.00001

void sc_fill_initial_approximation(double *approximation, size_t N);

void sc_fill_matrix_vector(double **matrix, double *vector, size_t N);

bool sc_is_solved(const double **matrix, const double *vector, double *curr_approximation, size_t N);

void sc_multiplication_matrix_vector(const double **matrix, const double *vector, double *res, size_t N);

void sc_subtracting_vectors(double *curr, const double *vector, size_t N);

void sc_multiplication_tau_vector(const double *vector, double *result, size_t N);

double sc_get_vector_sqrt(const double *vector, size_t N);

void sc_preparation_perfomance_free(size_t N);

void sc_solve_equations(const double **matrix, const double *vector, double *initial_approximation, size_t N);

void sc_print_result(double *res, size_t N);

#endif //SECOND_SCHEDULE_PROGRAM_H
