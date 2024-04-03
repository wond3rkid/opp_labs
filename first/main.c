#include <stdio.h>
#include <mpi/mpi.h>
#include <malloc.h>
#include <stdbool.h>

#define Epsilon 0.00001
#define Tau 0.00001

int rank, size;
int *chunk_array, *shift_array;
int N = 10;

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

void calculate_shift_array() {
    for (int i = 0; i < size; i++) {
        shift_array[i] = calculate_shift(i);
    }
}

void calculate_chunk_array() {
    for (int i = 0; i < size; i++) {
        chunk_array[i] = calculate_chunk_size(i);
    }
}

void fill_all_data(double *matrix, double *vector, double *approximation) {
    int curr_chunk = chunk_array[rank];
    int curr_shift = shift_array[rank];
    for (int i = 0; i < curr_chunk; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i * N + j] = (i == (j - curr_shift)) ? 2 : 1;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < N; i++) {
        approximation[i] = 0;
        vector[i] = N + 1;
    }
}


void parallel_print_result(double *vector) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (0 == rank) {
        for (int i = 0; i < N; i++) {
            printf("%1.1f ", vector[i]);
        }
    }
}


// return res or input in func?
double get_vector_sqrt(const double *vector) {
    double res;
    double chunk_res = 0;
    for (int i = 0; i < N; i++) {
        chunk_res += vector[i] * vector[i];
    }
    MPI_Allreduce(&chunk_res, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return res;
}

void mult_tau_vector(const double *vector, double *result) {
    for (int i = 0; i < N; i++) {
        result[i] = vector[i] * Tau;
    }
}

double *mult_matrix_vector(const double *matrix, const double *vector) {
    int curr_chunk = chunk_array[rank];
    double *result = calloc(curr_chunk, sizeof(double));
    for (int i = 0; i < curr_chunk; i++) {
        for (int j = 0; j < N; j++) {
            result[i] += matrix[i * N + j] * vector[j];
        }
    }
    return result;
}

double *subt_vectors(double *curr, const double *vector) {
    double *res = calloc(N, sizeof(double));
    for (int i = 0; i < N; i++) {
        res[i] = curr[i] - vector[i];
    }
    return res;
}


bool is_solved(double *den_vec, double *num_vec) {
    return get_vector_sqrt(den_vec) / get_vector_sqrt(num_vec) < Epsilon;
}

double *solve_equations(double *matrix, double *vector, double *approximation) {

    double *result = calloc(N, sizeof(double));
    bool is_done = is_solved(approximation, vector);
    fprintf(stderr, "1 IS SOLVED: %f \n", get_vector_sqrt(approximation) / get_vector_sqrt(vector));
    while (!is_done) {
        fprintf(stderr, " now is %d rank: %d \n", is_done, rank);
        // Ax Ax-b *tau x-tau(Ax-b)
        double *tmp_Ax = mult_matrix_vector(matrix, approximation);
        double *Ax = calloc(N, sizeof(double));
        MPI_Allgatherv(tmp_Ax, chunk_array[rank], MPI_DOUBLE, Ax, chunk_array, shift_array, MPI_DOUBLE, MPI_COMM_WORLD);

        double *Ax_b = subt_vectors(Ax, vector);
        fprintf(stderr, "IS SOLVED: %f \n", get_vector_sqrt(Ax_b) / get_vector_sqrt(vector));
        is_done = is_solved(Ax_b, vector);

        double *tAxminB = calloc(N, sizeof(double));
        mult_tau_vector(Ax_b, tAxminB);

        result = subt_vectors(result, tAxminB);
        printf("\n");
    }
    parallel_print_result(result);
    return result;
}


int main(int argc, char **argv) {
    double *matrix = malloc(sizeof(*matrix) * N * N);
    double *vector = malloc(sizeof(vector) * N);
    double *initial_approximation = malloc(sizeof(initial_approximation) * N);
    chunk_array = calloc(size, sizeof(int));
    shift_array = calloc(size, sizeof(int));
    double start, end;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    fill_all_data(matrix, vector, initial_approximation);
    calculate_chunk_array();
    calculate_shift_array();

    start = MPI_Wtime();
    solve_equations(matrix, vector, initial_approximation);
    end = MPI_Wtime();
    printf("Time taken for program: %f curr process: %d \n", end - start, rank);

    MPI_Finalize();

    free(matrix);
    free(vector);
    free(initial_approximation);
    free(chunk_array);
    free(shift_array);
    return MPI_SUCCESS;
}