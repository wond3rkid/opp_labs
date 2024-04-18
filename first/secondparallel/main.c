#include <stdio.h>
#include <malloc.h>
#include <mpi.h>
#include <math.h>

#define Epsilon 0.00001
#define Tau 0.00001
#define N 10
#define MATRIX_TAG 100
int rank, size;
int *shift_array, *chunk_array;

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
    int chunk_size = calculate_chunk_size(rank);
    double *vector = malloc(sizeof(double) * chunk_size);
    for (int i = 0; i < chunk_size; i++) {
        vector[i] = N + 1;
    }
    return vector;
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


void parallel_print_vector(const double *vector) {
    int chunk_size = calculate_chunk_size(rank);
    for (int proc = 0; proc < size; proc++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) {
            for (int i = 0; i < chunk_size; i++) {
                printf("%1.1f ", vector[i]);
            }
            printf("\n");
        }
    }
}

double get_vector_sqrt(const double *vector) {
    double chunk_res = 0;
    double res;
    for (int i = 0; i < chunk_array[rank]; i++) {
        chunk_res += vector[i] * vector[i];
    }
    MPI_Allreduce(&chunk_res, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return pow(res, 0.5);
}

double *mult_tau_vector(const double *vector) {
    double *res = malloc(chunk_array[rank] * sizeof(double));
    for (int i = 0; i < chunk_array[rank]; i++) {
        res[i] = vector[i] * Tau;
    }
    return res;
}

//double* multiply_matrix_by_vector(double* matrix, double* vector) {
//    int vector_length = c_recvcounts[rank];
//    int offset = c_displs[rank];
//    int incoming_process_data = 0;
//
//    double* result = (double*) calloc(vector_length, sizeof(double));
//    double* v = (double*) calloc(N/size + 1, sizeof(double));
//    for (int i = 0; i < vector_length; i++) {
//        v[i] = vector[i];
//    }
//    for (int process = 0; process < size; process++) {
//        // index of current part of vector
//        incoming_process_data = (rank + process) % size;
//
//        for (int i = 0; i < c_recvcounts[rank]; i++) {
//            for (int j = 0; j < c_recvcounts[incoming_process_data]; j++) {
//                result[i] += matrix[i * N + j + c_displs[incoming_process_data]] * v[j];
//            }
//        }
//        // switch vector, (vector length, vector offset) between processes
//        MPI_Sendrecv_replace(v, N/size + 1, MPI_DOUBLE, (rank+1) % size, TAG_SEND_MATRIX,
//                             (rank-1) % size, TAG_SEND_MATRIX, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//    }
//    return result;
//}

double *mult_matrix_vector(const double *matrix, const double *vector) {
    int chunk_size = chunk_array[rank];
    int tmp_proc;

    double *res = calloc(chunk_size, sizeof(double));
    // ? chunk_size
    double *tmp = malloc(chunk_size * sizeof(double));

    for (int i = 0; i < chunk_size; i++) {
        tmp[i] = vector[i];
    }
    for (int proc = 0; proc < size; proc++) {
        tmp_proc = (rank + proc) % N;

        for (int i = 0; i < shift_array[rank]; i++) {
            for (int j = 0; j < shift_array[tmp_proc]; j++) {
                res[i] += matrix[i * N + j + shift_array[tmp_proc]] * tmp[j];
            }
        }
        MPI_Sendrecv_replace(tmp, chunk_size, MPI_DOUBLE, (rank + 1) % size, MATRIX_TAG, (rank) % size,
                             MATRIX_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }
    return res;
}

double *subt_vectors(const double *curr, const double *vector) {
    double *res = calloc(chunk_array[rank], sizeof(double));
    for (int i = 0; i < chunk_array[rank]; i++) {
        res[i] = curr[i] - vector[i];
    }
    return res;
}

bool is_solved(double *axb, double *b) {
    double axb_sqrt = get_vector_sqrt(axb);
    double b_sqrt = get_vector_sqrt(b);
    return axb_sqrt / b_sqrt < Epsilon;
}

void solve_equations(double *matrix, double *vector) {
    double *result = calloc(N, sizeof(double));

    while (true) {
        double *Ax = mult_matrix_vector(matrix, result);
        double *Axb = subt_vectors(Ax, vector);

        if (is_solved(Axb, vector)) {
            printf("System of linear algebraic equations was solved. Check the result:\n");
            parallel_print_vector(result);
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
    chunk_array = calculate_chunk_array();
    shift_array = calculate_shift_array();
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        printf("Input matrix : \n");
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