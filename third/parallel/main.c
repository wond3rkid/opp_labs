


#include <mpi.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>

#define N1 3600
#define N2 1200
#define N3 2400

void create_known_filling(double *matrix, int r, int c) {
    assert(matrix != NULL);
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            matrix[i * r + c] = i == j ? 1 : 0;
        }
    }
}

void create_rand_filling(double *matrix, int r, int c) {
    assert(matrix != NULL);
    srand(time(NULL));
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            matrix[i * r + c] = rand() % 5;
        }
    }
}

void multiplicate(const double *first, const double *second, double *result, int a, int b, int c) {
    assert(first != NULL);
    assert(second != NULL);
    assert(result != NULL);
    for (int i = 0; i < a; i++) {
        for (int j = 0; j < c; j++) {
            result[i * c + j] = 0;
            for (int k = 0; k < b; ++k) {
                result[i * c + j] += first[i * b + k] * second[k * c + j];
            }
        }
    }
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "Input error");
        return EXIT_FAILURE;
    }
    // инициализация mpi
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm comm, row_comm, column_comm;
    //
    int x, y;
    x = atoi(argv[1]);
    y = atoi(argv[2]);
    fprintf(stdout, "Your input x - %d and y - %d \n", x, y);
    int dims[] = {x, y};
    int periods[] = {0, 0};
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &comm);

    int remain_dims[2];
    remain_dims[0] = 0;
    remain_dims[1] = 1;
    MPI_Cart_sub(comm, remain_dims, &row_comm);

    remain_dims[0] = 1;
    remain_dims[1] = 0;
    MPI_Cart_sub(comm, remain_dims, &column_comm);

    int column_size, row_size;
    MPI_Comm_size(column_comm, &column_size);
    MPI_Comm_size(row_comm, &row_size);

    int coords[2];
    MPI_Cart_coords(comm, rank, 2, coords);
    double *a_part, *b_part, *c_part;
    a_part = (double *) calloc(N1 * N2 / x, sizeof(double));
    b_part = (double *) calloc(N2 * N3 / y, sizeof(double));
    c_part = (double *) calloc(N1 * N2 / x / y, sizeof(double));

    double *a, *b, *c;
    double start_time = MPI_Wtime();
    if (rank == 0) {
        a = (double *) calloc(sizeof(double), N1 * N2);
        b = (double *) calloc(sizeof(double), N2 * N3);
        c = (double *) calloc(sizeof(double), N1 * N3);
        if (a == NULL || b == NULL || c == NULL) {
            fprintf(stderr, "Variables a, b, or c are not initialized\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        create_known_filling(a, N1, N2);
        create_known_filling(b, N2, N3);
    }
    if (coords[1] == 0) {
        MPI_Scatter(a, N1 * N2 / x, MPI_DOUBLE, a_part, N1 * N2 / x, MPI_DOUBLE, 0, column_comm);
    }

    MPI_Datatype col_comm;
    MPI_Type_vector(N2, N3 / y, N3, MPI_DOUBLE, &col_comm);
    MPI_Type_commit(&col_comm);
    MPI_Datatype c_type;

    MPI_Type_vector(N1 / x, N3 / y, N3, MPI_DOUBLE, &c_type);
    MPI_Type_commit(&c_type);

    if (rank == 0) {
        for (int i = 1; i < row_size; i++) {
            MPI_Send(b + i * N3 / y, 1, col_comm, i, 0, row_comm);
        }

        for (int i = 0; i < N2; i++) {
            for (int j = 0; j < N3 / y; j++) {
                b_part[i * N3 / y + j] = b[i * N3 + j];
            }
        }
    } else if (coords[0] == 0) {
        MPI_Recv(b_part, N3 * N2 / y, MPI_DOUBLE, 0, 0, row_comm, MPI_STATUS_IGNORE);
    }

    MPI_Bcast(a_part, N1 * N2 / x, MPI_DOUBLE, 0, row_comm);
    MPI_Bcast(b_part, N2 * N3 / y, MPI_DOUBLE, 0, column_comm);
    multiplicate(a_part, b_part, c_part, N1 / x, N3 / y, N2);

    if (rank != 0) {
        MPI_Send(c_part, N3 * N1 / y / x, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    } else {
        for (int i = 0; i < N1 / x; i++) {
            for (int j = 0; j < N3 / y; j++) {
                c[i * N3 + j] = c_part[i * N3 / y + j];
            }
        }
        for (int i = 0; i < x; i++) {
            for (int j = 0; j < y; j++) {
                if (i != 0 || j != 0)
                    MPI_Recv(c + i * N3 * N1 / x + j * N3 / y, 1, c_type, i * y + j, 1, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
            }
        }
        //free(a);
        //free(b);
        //free(c);
    }
    double finish_time = MPI_Wtime();
    fprintf(stdout, "proc %d time %lf \n ", rank, finish_time - start_time);
    //free(a_part);
    //free(b_part);
    //free(c_part);
    MPI_Type_free(&col_comm);
    MPI_Type_free(&c_type);
    MPI_Finalize();
    return EXIT_SUCCESS;
}

