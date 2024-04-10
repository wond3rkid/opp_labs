#include <stdio.h>
#include <malloc.h>
#include <mpi.h>
#include <math.h>

#define Epsilon 0.00001
#define Tau 0.00001
#define N 10000
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

double *create_matrix() {
    double *matrix = malloc(sizeof(double) * N * N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i * N + j] = (i == j) ? 2 : 1;
        }
    }
    return matrix;
}

double *create_vector() {
    double *vector = malloc(sizeof(double) * N);
    for (int i = 0; i < N; i++) {
        vector[i] = N + 1;
    }
    return vector;
}

void parallel_print_result(const double *vector) {
    if (rank == 0) {
        printf("System of linear algebraic equations was solved. Check the result:\n");
        for (int i = 0; i < N; i++) {
            printf("%1.1f ", vector[i]);
        }
        printf("\n");
    }
}

double get_vector_sqrt(const double *vector) {
    double res = 0;
    int chunk_size = calculate_chunk_size(rank);
    for (int i = 0; i < chunk_size; i++) {
        res += vector[i] * vector[i];
    }
    return pow(res, 0.5);
}

double *mult_tau_vector(const double *vector) {
    int chunk_size = calculate_chunk_size(rank);
    double *res = malloc(chunk_size * chunk_size);
    for (int i = 0; i < chunk_size; i++) {
        res[i] = vector[i] * Tau;
    }
    return res;
}

double *mult_matrix_vector(const double *matrix, const double *vector) {
    int chunk_size = calculate_chunk_size(rank);
    double *result = calloc(chunk_size, sizeof(double));
    for (int i = 0; i < chunk_size; i++) {
        for (int j = 0; j < N; j++) {
            result[i] += matrix[i * N + j] * vector[j];
        }
    }
    return result;
}

double *subt_vectors(const double *curr, const double *vector) {
    int chunk_size = calculate_chunk_size(rank);
    double *res = calloc(chunk_size, sizeof(double));
    for (int i = 0; i < chunk_size; i++) {
        res[i] = curr[i] - vector[i];
    }
    return res;
}

bool is_solved(double *axb, double *b) {
    double axb_sqrt_part = get_vector_sqrt(axb);
    double axb_sqrt = 0;
    MPI_Allreduce(&axb_sqrt_part, &axb_sqrt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double b_sqrt_part = get_vector_sqrt(b);
    double b_sqrt = 0;
    MPI_Allreduce(&b_sqrt_part, &b_sqrt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return axb_sqrt / b_sqrt < Epsilon;
}

void solve_equations(double *matrix, double *vector) {
    double *result = calloc(N, sizeof(double));
    int *chunk_array = calculate_chunk_array();
    int *shift_array = calculate_shift_array();
    while (true) {
        double *Ax_tmp = mult_matrix_vector(matrix, result);
        double *Ax = calloc(N, sizeof(double));
        MPI_Allgatherv(Ax_tmp, chunk_array[rank], MPI_DOUBLE, Ax, chunk_array, shift_array, MPI_DOUBLE, MPI_COMM_WORLD);

        double *Axb_tmp = subt_vectors(Ax, vector);
        double *Axb = calloc(N, sizeof(double));
        MPI_Allgatherv(Axb_tmp, chunk_array[rank], MPI_DOUBLE, Axb, chunk_array, shift_array, MPI_DOUBLE, MPI_COMM_WORLD);

        if (is_solved(Axb, vector)) {
            parallel_print_result(result);
            free(Ax_tmp);
            free(Ax);
            free(Axb);
            free(Axb_tmp);
            break;
        }
        double *tAxb_tmp = mult_tau_vector(Axb);
        double *tAxb = calloc(sizeof(double), N);
        MPI_Allgatherv(tAxb_tmp, chunk_array[rank], MPI_DOUBLE, tAxb, chunk_array, shift_array, MPI_DOUBLE,
                       MPI_COMM_WORLD);

        double *result_tmp = subt_vectors(result, tAxb);

        MPI_Allgatherv(result_tmp, chunk_array[rank], MPI_DOUBLE, result, chunk_array, shift_array, MPI_DOUBLE, MPI_COMM_WORLD);
        free(tAxb);
        free(tAxb_tmp);
        free(result_tmp);
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