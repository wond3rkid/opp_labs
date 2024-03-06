#ifndef SECOND_PARALLEL_PROGRAM_H
#define SECOND_PARALLEL_PROGRAM_H
#include <omp.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdbool.h>

#define Epsilon 0.00001
#define Tau 0.00001

void p_fill_initial_approximation(double *approximation, size_t N);

void p_fill_matrix_vector(double **matrix, double *vector, size_t N);

bool p_is_solved(const double **matrix, const double *vector, double *curr_approximation, size_t N);

double *p_multiplication_matrix_vector(const double **matrix, const double *vector, double *res, size_t N);

void p_subtracting_vectors(double *curr, const double *vector, size_t N);

double *p_multiplication_tau_vector(const double *vector, size_t N);

double p_get_vector_sqrt(const double *vector, size_t N);

double *p_get_next_x(const double **matrix, const double *vector, double *curr_approximation, size_t N);

void p_preparation_perfomance_free(size_t N);

void p_solve_equations(const double **matrix, const double *vector, double *initial_approximation, size_t N);

void p_print_result(double *res, size_t N);

#endif //SECOND_PARALLEL_PROGRAM_H
