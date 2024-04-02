#include "program.h"

void fill_all_data(double *matrix, double *vector, double *approximation) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i * N + j] = (i == j) ? 2 : 1;
        }
    }
    for (int i = 0; i < N; i++) {
        approximation[i] = 0;
        vector[i] = N + 1;
    }
}