#include <stdio.h>
#include <mpi/mpi.h>
#include <malloc.h>
#include "program.h"

int main(int argc, char **argv) {
    int
            N_ = 1000;
    double *matrix = malloc(sizeof(*matrix) * N_ * N_);
    double *vector = malloc(sizeof(vector) * N_);
    double *initial_approximation = malloc(sizeof(initial_approximation) * N_);
    MPI_Init(&argc, &argv);
    int size, rank;
    double start = 0, end;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        fill_all_data(matrix, vector, initial_approximation);
        start = MPI_Wtime();
    }


    end = MPI_Wtime();
    printf("Time taken for program: %f \n", end - start);
    MPI_Finalize();
}