#include <stdio.h>
#include <mpi/mpi.h>
#include <malloc.h>

#define Epsilon 0.00001
#define Tau 0.00001
#define N 10
int rank, size;
int *chunk_array, *shift_array;

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

void parallel_print_result(double *vector) {
    MPI_Barrier(MPI_COMM_WORLD);
//    if (0 == rank) {
    for (int i = 0; i < N; i++) {
        printf("%1.1f ", vector[i]);
    }
//    }
}

double get_vector_sqrt(const double *vector) {
    double res = 0;
    for (int i = 0; i < N; i++) {
        res += vector[i] * vector[i];
    }
    return res;
}

void mult_tau_vector(const double *vector, double *res) {
    for (int i = 0; i < N; i++) {
        res[i] = vector[i] * Tau;
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

double *subt_vectors(const double *curr, const double *vector) {
    double *res = calloc(N, sizeof(double));
    for (int i = 0; i < N; i++) {
        res[i] = curr[i] - vector[i];
    }
    return res;
}

void subt_in_place(double *result, const double *minus) {

    for (int i = 0; i < N; i++) {
        result[i] -= minus[i];
    }
}

bool another_is_solved(const double *matrix, const double *vector, const double *curr_approximation) {
    double numerator_sqrt = get_vector_sqrt(vector);
    double *denominator = mult_matrix_vector(matrix, curr_approximation);
    double *full_den = calloc(N, sizeof(double));
    MPI_Allgatherv(denominator, chunk_array[rank], MPI_DOUBLE, full_den, chunk_array, shift_array, MPI_DOUBLE,
                   MPI_COMM_WORLD);
    subt_in_place(full_den, vector);
    double denominator_sqrt = get_vector_sqrt(full_den);
    free(denominator);
    free(full_den);
    return denominator_sqrt / numerator_sqrt < Epsilon;
}


void solve_equations(double *matrix, double *vector, double *approximation) {
    do {
        double *tmp;
        tmp = mult_matrix_vector(matrix, approximation);
        double *Ax = calloc(N, sizeof(double));
        MPI_Allgatherv(tmp, chunk_array[rank], MPI_DOUBLE, Ax, chunk_array, shift_array, MPI_DOUBLE, MPI_COMM_WORLD);
        subt_in_place(tmp, vector);
        double *curr = calloc(N, sizeof(double));
        mult_tau_vector(tmp, curr);
        subt_in_place(approximation, curr);
        parallel_print_result(approximation);
        free(curr);
        free(Ax);
        free(tmp);

    } while (!another_is_solved(matrix, vector, approximation));
    parallel_print_result(approximation);
}

int main(int argc, char **argv) {
    double *matrix = malloc(sizeof(*matrix) * N * N);
    double *vector = malloc(sizeof(vector) * N);
    double *initial_approximation = malloc(sizeof(initial_approximation) * N);
    double start, end;
    fill_all_data(matrix, vector, initial_approximation);
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Barrier(MPI_COMM_WORLD);
//    if (rank == 0) { //FIXME
    chunk_array = calloc(size, sizeof(int));
    shift_array = calloc(size, sizeof(int));
    calculate_chunk_array();
    calculate_shift_array();
//    }
    fprintf(stderr, "sizeof size: %d rank : %d curr chunk: %d \n", size, rank, chunk_array[rank]);
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