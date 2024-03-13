#ifndef SECOND_PARALLEL_SECTION_PROGRAM_H
#define SECOND_PARALLEL_SECTION_PROGRAM_H

#include <omp.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define Epsilon 0.00001
#define Tau 0.00001

void s_fill_initial_approximation(double *approximation, size_t N);

void s_fill_matrix_vector(double **matrix, double *vector, size_t N);

bool s_is_solved(const double **matrix, const double *vector, double *curr_approximation, size_t N);

void s_multiplication_matrix_vector(const double **matrix, const double *vector, double *res, size_t N);

void s_subtracting_vectors(double *curr, const double *vector, size_t N);

void s_multiplication_tau_vector(const double *vector, double *result, size_t N);

double s_get_vector_sqrt(const double *vector, size_t N);

void s_get_next_x(const double **matrix, const double *vector, double *curr_approximation, size_t N);

void s_preparation_perfomance_free(size_t N);

void s_solve_equations(const double **matrix, const double *vector, double *initial_approximation, size_t N);

void p_print_result(double *res, size_t N);

#endif //SECOND_PARALLEL_SECTION_PROGRAM_H
