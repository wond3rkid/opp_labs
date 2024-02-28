#ifndef SECOND_PROCESS_H
#define SECOND_PROCESS_H

#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdbool.h>

#define Epsilon 0,00001
#define Tau 0.01

// заполнение вектора с начальным приближением
void fill_vector_initial_approximation(double *approximation, size_t N);

// заполнение матрицы А и вектора b по условию задания
void fill_matrix_vector(double **matrix, double *vector, size_t N);

// проверка решена ли система == достигнуто ли искомое приближение эпсилон
bool is_solved(double **matrix, double *vector, double *curr_approximation, size_t N);

// перемножение матрицы на вектор и возвращение нового вектора

double *multiplication_matrix_vector(double **matrix, double *curr_approximation, size_t N);

// вычитание из вектора с результатом перемножения матрицы и вектора вектора свободных значений b
void subtracting_vectors(double *curr, double *vector, size_t N);

double  **multiplication_tau_matrix(double **matrix, size_t N;

double get_vector_sqrt(double *vector, size_t N);

double *iterationMethod(double **matrix_A, double *matrix_B, size_t N);

void do_algorithm(size_t N);

void print_result(double *res, size_t N);

#endif //SECOND_PROCESS_H
