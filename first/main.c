#include <stdio.h>
#include <mpi/mpi.h>
#include <malloc.h>

int rank, size;
int N = 100;

int calculate_chunk_size(int rank, int size) {
    int rest = N % size;
    int add = (rank < rest) ? 1 : 0;
    return (int) (N / size) + add;
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

int calculate_shift(int rank, int size) {
    int shift = 0;
    for (int i = 0; i < rank; i++) {
        shift += calculate_chunk_size(i, size);
    }
    return shift;
}

void calculate_shift_array(int *shift_array, int size) {
    for (int i = 0; i < size; i++) {
        shift_array[i] = calculate_shift(i, size);
    }
}

void calculate_chunk_array(int *chunk_array, int size) {
    for (int i = 0; i < size; i++) {
        chunk_array[i] = calculate_chunk_size(i, size);
    }
}

void parallel_print_result(double *vector, int size, int *chunk_sizes) {
    for (int proc = 0; proc < size; proc++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (proc == rank) {
            for (int i = 0; i < chunk_sizes[rank]; i++) {
                printf("%1.1f", vector[i]);
            }
        }
    }
}


int main(int argc, char **argv) {
    int array_size = 1000;
    double *matrix = malloc(sizeof(*matrix) * array_size * array_size);
    double *vector = malloc(sizeof(vector) * array_size);
    double *initial_approximation = malloc(sizeof(initial_approximation) * array_size);
    int *chunk_array, *shift_array;

    MPI_Init(&argc, &argv);
    double start = 0, end;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("rank : %d size: %d \n", rank, size);

    if (rank == 0) {
        //preparation
        chunk_array = calloc(sizeof(int), size);
        shift_array = calloc(sizeof(int), size);
        fill_all_data(matrix, vector, initial_approximation);
        calculate_chunk_array(chunk_array, size);
        calculate_shift_array(shift_array, size);
        start = MPI_Wtime();
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        end = MPI_Wtime();
        printf("Time taken for program: %f \n", end - start);
    }

    MPI_Finalize();
    return MPI_SUCCESS;
}