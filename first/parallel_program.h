#ifndef FIRST_PARALLEL_PROGRAM_H
#define FIRST_PARALLEL_PROGRAM_H

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdbool.h>

#define Epsilon 0.00001
#define Tau 0.00001

// заполнение вектора с начальным приближением
void fill_vector_initial_approximation(double *approximation, size_t N);

// заполнение матрицы А и вектора b по условию задания
void fill_matrix_vector(double **matrix, double *vector, size_t N);

// проверка решена ли система == достигнуто ли искомое приближение эпсилон
bool is_solved(const double **matrix, const double *vector, double *curr_approximation, size_t N);

// перемножение матрицы на вектор и возвращение нового вектора
void multiplication_matrix_vector(const double **matrix, const double *vector, double *res, size_t N);

// вычитание из вектора с результатом перемножения матрицы и вектора вектора свободных значений b
void subtracting_vectors(double *curr, const double *vector, size_t N);

// умножение матрицы на число Тау. Возвращает новую матрицу
void multiplication_tau_vector(const double *vector, double *result, size_t N);

// считает норму вектора
double get_vector_sqrt(const double *vector, size_t N);

void do_algorithm(size_t N);

void solve_equations(const double **matrix, const double *vector, double *initial_approximation, size_t N);

void print_result(double *res, size_t N);

#endif //FIRST_PARALLEL_PROGRAM_H
