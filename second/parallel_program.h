#ifndef SECOND_PARALLEL_PROGRAM_H
#define SECOND_PARALLEL_PROGRAM_H
#include <omp.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdbool.h>

#define Epsilon 0.00001
#define Tau 0.00001

void parallel_fill_initial_approximation(double *approximation, size_t N);

void parallel_fill_matrix_vector(double **matrix, double *vector, size_t N);

bool is_solved_parallel(const double **matrix, const double *vector, double *curr_approximation, size_t N);

double *parallel_multiplication_matrix_vector(const double **matrix, double *curr_approximation, size_t N);

void parallel_subtracting_vectors(double *curr, const double *vector, size_t N);

double *parallel_multiplication_tau_vector(const double *vector, size_t N);

double get_vector_sqrt_parallel(const double *vector, size_t N);

double *get_next_x_parallel(const double **matrix, const double *vector, double *curr_approximation, size_t N);

void preparation_perfomance_free_parallel(size_t N);

void parallel_solve_equations(const double **matrix, const double *vector, double *initial_approximation, size_t N);

void print_result_for_parallel(double *res, size_t N);

#endif //SECOND_PARALLEL_PROGRAM_H
