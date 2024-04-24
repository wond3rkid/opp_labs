#include <stdio.h>
#include <malloc.h>
#include <mpi.h>
#include <math.h>
#include <unistd.h>

#define Epsilon 0.00001
#define Tau 0.00001
#define N 10
int rank, size;

int calculate_chunk_size(int curr_rank) {
    int rest = N % size;
    int add = (curr_rank < rest) ? 1 : 0;
    return (int) (N / size) + add;
}

int calculate_shift(int curr_rank) {
    int shift = 0;
    for (int i = 0; i < curr_rank; i++) {
        shift += calculate_chunk_size(i);
    }
    return shift;
}

int *calculate_shift_array() {
    int *curr_shift_array = malloc(sizeof(int) * size);
    for (int i = 0; i < size; i++) {
        curr_shift_array[i] = calculate_shift(i);
    }
    return curr_shift_array;
}

int *calculate_chunk_array() {
    int *curr_chunk_array = malloc(sizeof(int) * size);
    for (int i = 0; i < size; i++) {
        curr_chunk_array[i] = calculate_chunk_size(i);
    }
    return curr_chunk_array;
}

void parallel_print_matrix(const double *matrix) {
    int chunk_size = calculate_chunk_size(rank);
    for (int proc = 0; proc < size; proc++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (proc == rank) {
            for (int i = 0; i < chunk_size; i++) {
                for (int j = 0; j < N; j++) {
                    printf("%0.0f ", matrix[i * N + j]);
                }
                printf("\n");
            }
        }
    }
}

double *create_matrix() {
    int chunk_size = calculate_chunk_size(rank);
    int shift = calculate_shift(rank);
    double *matrix_chunk = malloc(sizeof(double) * N * chunk_size);
    for (int i = 0; i < chunk_size; i++) {
        for (int j = 0; j < N; j++) {
            matrix_chunk[i * N + j] = (i == (j - shift)) ? 2 : 1;
        }
    }
    return matrix_chunk;
}

double *create_vector() {
    double *vector = malloc(sizeof(double) * N);
    for (int i = 0; i < N; i++) {
        vector[i] = N + 1;
    }
    return vector;
}

void parallel_print_vector(const double *vector) {
    if (rank == 0) {
        for (int i = 0; i < N; i++) {
            printf("%1.1f ", vector[i]);
        }
        printf("\n");
    }
}


double get_vector_sqrt(const double *vector) {
    double res = 0;
    for (int i = 0; i < N; i++) {
        res += vector[i] * vector[i];
    }
    return pow(res, 0.5);
}

double *mult_tau_vector(const double *vector) {
    double *res = malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        res[i] = vector[i] * Tau;
    }
    return res;
}

double *mult_matrix_vector(const double *matrix, const double *vector) {
    int curr_chunk = calculate_chunk_size(rank);
    double *result = calloc(curr_chunk, sizeof(double));
    for (int i = 0; i < curr_chunk; i++) {
        for (int j = 0; j < N; j++) {
            result[i] += matrix[i * N + j] * vector[j];
        }
    }
    return result;
}

double *subt_vectors(const double *curr, const double *vector) {
    double *res = calloc(N, sizeof(double));
    for (int i = 0; i < N; i++) {
        res[i] = curr[i] - vector[i];
    }
    return res;
}


bool is_solved(double *axb, double *b) {
    return get_vector_sqrt(axb) / get_vector_sqrt(b) < Epsilon;
}

void solve_equations(double *matrix, double *vector) {
    double *result = calloc(N, sizeof(double));
    int *chunk_array = calculate_chunk_array();
    int *shift_array = calculate_shift_array();
    while (true) {
        double *tmp = mult_matrix_vector(matrix, result);
        double *Ax = calloc(N, sizeof(double));
        MPI_Allgatherv(tmp, chunk_array[rank], MPI_DOUBLE, Ax, chunk_array, shift_array, MPI_DOUBLE, MPI_COMM_WORLD);
        double *Axb = subt_vectors(Ax, vector);
        if (is_solved(Axb, vector)) {
            if (rank == 0) {
                printf("Result of solving equations: \n");
            }
            parallel_print_vector(result);
            free(tmp);
            free(Ax);
            free(Axb);
            break;
        }
        double *tAxb = mult_tau_vector(Axb);
        result = subt_vectors(result, tAxb);
        free(tAxb);
    }
    free(result);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Barrier(MPI_COMM_WORLD);

    double *vector = create_vector();
    double *matrix = create_matrix();
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        printf("Input matrix N = %d: \n", N);
    }
    parallel_print_matrix(matrix);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        printf("Input vector: \n");
    }
    parallel_print_vector(vector);
    double start, end;

    start = MPI_Wtime();
    solve_equations(matrix, vector);
    end = MPI_Wtime();

    printf("Time taken for %d process : %f sec\n", rank, end - start);
    free(matrix);
    free(vector);

    MPI_Finalize();
    return MPI_SUCCESS;
}