#include <stdio.h>
#include <malloc.h>
#include <mpi.h>
#include <stdbool.h>

#define Epsilon 0.00001
#define Tau  0.00001
#define N 25000
#define MATRIX_TAG 100

int rank, size;
int *displs;
int *recvcounts;

int calculate_chunk(int curr_rank) {
    int part = N / size;
    int rest = N % size;
    return part + (curr_rank < rest ? 1 : 0);
}

int calculate_shift(int curr_rank) {
    int shift = 0;
    for (int i = 0; i < curr_rank; i++) {
        shift += calculate_chunk(i);
    }
    return shift;
}

int *create_displs() {
    int *temp = (int *) calloc(size, sizeof(int));
    for (int i = 0; i < size; i++) {
        temp[i] = calculate_shift(i);
    }
    return temp;
}

int *create_revcounts() {
    int *temp = (int *) calloc(size, sizeof(int));
    for (int i = 0; i < size; i++) {
        temp[i] = calculate_chunk(i);
    }
    return temp;
}

double *create_matrix() {
    int chunk = recvcounts[rank];
    int shift = displs[rank];
    double *matrix_chunk = (double *) calloc(N * chunk, sizeof(double));
    for (int i = 0; i < chunk; i++) {
        for (int j = 0; j < N; j++) {
            matrix_chunk[i * N + j] = (i == (j - shift)) ? 2 : 1;
        }
    }
    return matrix_chunk;
}

double *create_vector() {
    int chunk = recvcounts[rank];
    double *vector_chunk = (double *) calloc(chunk, sizeof(double));
    for (int i = 0; i < chunk; i++) {
        vector_chunk[i] = N + 1;
    }
    return vector_chunk;
}

void parallel_print_matrix(const double *matrix) {
    int chunk = calculate_chunk(rank);
    for (int proc = 0; proc < size; proc++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (proc == rank) {
            for (int i = 0; i < chunk; i++) {
                for (int j = 0; j < N; j++) {
                    printf("%0.0f ", matrix[i * N + j]);
                }
                printf("\n");
            }
        }
    }
}


void parallel_print_vector(const double *vector) {
    int chunk = calculate_chunk(rank);
    for (int proc = 0; proc < size; proc++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == proc) {
            for (int i = 0; i < chunk; i++) {
                printf("%1.1f ", vector[i]);
            }
            printf("\n");
        }
    }
}

double get_vector_sqrt(const double *vector) {
    double chunk_res = 0;
    double res;
    for (int i = 0; i < recvcounts[rank]; i++) {
        chunk_res += vector[i] * vector[i];
    }
    MPI_Allreduce(&chunk_res, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return res;
}

double *subtraction(const double *minuend, const double *subtrahend) {
    double *res = (double *) calloc(recvcounts[rank], sizeof(double));
    for (int i = 0; i < recvcounts[rank]; i++) {
        res[i] = minuend[i] - subtrahend[i];
    }
    return res;
}

double *multiplication_tau(const double *vector) {
    double *res = (double *) calloc(recvcounts[rank], sizeof(double));
    for (int i = 0; i < recvcounts[rank]; i++) {
        res[i] = vector[i] * Tau;
    }
    return res;
}

double *multiplication_matrix_vector(const double *matrix, const double *vector) {
    int chunk = recvcounts[rank];
    int temp_size;

    double *tmp = (double *) calloc(N / size + 1, sizeof(double));
    double *res = (double *) calloc(chunk, sizeof(double));
    for (int i = 0; i < chunk; i++) {
        tmp[i] = vector[i];
    }
    for (int proc = 0; proc < size; proc++) {
        temp_size = (rank + proc) % size;
        for (int i = 0; i < recvcounts[rank]; i++) {
            for (int j = 0; j < recvcounts[temp_size]; j++) {
                res[i] += matrix[i * N + j + displs[temp_size]] * tmp[i];
            }
        }
        MPI_Sendrecv_replace(tmp, N / size + 1, MPI_DOUBLE, (rank + 1) % size, MATRIX_TAG,
                             (rank + size - 1) % size, MATRIX_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    return res;
}

void solve_equations(double *matrix, double *vector) {
    double b_sqrt = get_vector_sqrt(vector);
    double *result = (double *) calloc(N, sizeof(double));
    while (true) {
        double *Ax = multiplication_matrix_vector(matrix, result);
        double *Axb = subtraction(Ax, vector);
        double axb_sqrt = get_vector_sqrt(Axb);
        if (axb_sqrt / b_sqrt < Epsilon * Epsilon) {
            printf("System of linear algebraic equations was solved. Check the result for proc %d:\n", rank);
            //parallel_print_vector(result);
            free(Ax);
            free(Axb);
            break;
        }
        double *tAxb = multiplication_tau(Axb);
        result = subtraction(result, tAxb);
        free(tAxb);
    }
    free(result);
}


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    recvcounts = create_revcounts();
    displs = create_displs();

    double *vector = create_vector();
    double *matrix = create_matrix();
    double start = MPI_Wtime();
    solve_equations(matrix, vector);
    double end = MPI_Wtime();
    printf("Time taken for %d process : %f sec\n", rank, end - start);

    free(recvcounts);
    free(displs);
    free(matrix);
    free(vector);
    MPI_Finalize();
    return MPI_SUCCESS;
}