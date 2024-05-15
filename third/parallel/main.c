#include <mpi/mpi.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>

#define N1 3600
#define N2 3600
#define N3 3600

void create_known_filling(double *matrix, int r, int c) {
    assert(matrix != NULL);
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            matrix[i * r + c] = i == j ? 1 : 0;
        }
    }
}

double create_rand_filling(double *matrix, int r, int c) {
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
    // инициализация mpi,
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm comm, row_comm, column_comm;
    // инициализируем решетку размером 2
    int x, y;
    x = atoi(argv[1]);
    y = atoi(argv[2]);
    fprintf(stdout, "Your input x - %d and y - %d \n", x, y);
    int dims[] = {x, y};
    int periods[] = {0, 0};
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &comm);
    // создаем типы коммуникаторов для раздачи по строкам и столбцам. remain_dims -
    int remain_dims[2];
    remain_dims[0] = 0;
    remain_dims[1] = 1;
    MPI_Cart_sub(comm, remain_dims, &row_comm);
    // комм - решетка, разделение на две решетки. мы будто бы создаем роу_комм, но их на самом деле столько же сколько строк
    // разделение на строчные и колоночные коммуникаторы
    // под х 0 не учитываем это измерение по оси, у
    remain_dims[0] = 1;
    remain_dims[1] = 0;
    MPI_Cart_sub(comm, remain_dims, &column_comm);
    // сейчас у каждого столбца и строчки есть коммуникаторы. берем кол-во процессов в коммуникаторах
    int column_size, row_size;
    MPI_Comm_size(column_comm, &column_size);
    MPI_Comm_size(row_comm, &row_size);
    // координаты из решетки
    int coords[2];
    MPI_Cart_coords(comm, rank, 2, coords);
    // выделения частей для хранения матрицы
    double *a_part, *b_part, *c_part;
    a_part = calloc(N1 * N2 / x, sizeof(double));
    b_part = calloc(N2 * N3 / y, sizeof(double));
    c_part = calloc(N1 * N2 / x / y, sizeof(double));

    double *a, *b, *c;
    double start_time;
    // на нулевом процессе заполняем
    if (rank == 0) {
        a = (double *) calloc(sizeof(double), N1 * N2);
        b = (double *) calloc(sizeof(double), N2 * N3);
        c = (double *) calloc(sizeof(double), N1 * N3);
        create_known_filling(a, N1, N2);
        create_known_filling(b, N2, N3);
        start_time = MPI_Wtime();
    }

    if (coords[1] == 0) {
        // первый процесс своим коммуникатором всем последующим
        MPI_Scatter(a, N1 * N2 / x, MPI_DOUBLE, a_part, N1 * N2 / x, MPI_DOUBLE, 0, column_comm); //!!!!
    }
    /* Матрица А без проблем разрезается по полосам, с матрицей В проблемы, создаем свой дататайп*/
    MPI_Datatype col_comm;
    MPI_Type_vector(N2, N3 / y, N3, MPI_DOUBLE, &col_comm);
    MPI_Type_commit(&col_comm);
    MPI_Datatype c_type; // сбор данных с матрицы С в 1 процесс (миноры матрицы С)

    MPI_Type_vector(N1 / x, N3 / y, N3, MPI_DOUBLE, &c_type);
    MPI_Type_commit(&c_type);

    if (rank == 0) {
        for (int i = 1; i < row_size; i++) {
            // разделяем матрицу B по процессам
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
    // bcast делим парты по процессам
    MPI_Bcast(a_part, N1 * N2 / x, MPI_DOUBLE, 0, row_comm);
    MPI_Bcast(b_part, N2 * N3 / y, MPI_DOUBLE, 0, column_comm);
    // перемножаем миноры и собираем в 0 процесс, в сам себя не бросаем - это долго
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
                    // нулевой процесс не ждет отправки сам себе, c_part -> c
                    MPI_Recv(c + i * N3 * N1 / x + j * N3 / y, 1, c_type, i * y + j, 1, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
            }
        }
        double finish_time = MPI_Wtime();
        printf("%lf\n", finish_time - start_time + 0.000000001 * (finish_time - start_time));
        free(a);
        free(b);
        free(c);
    }
    free(a_part);
    free(b_part);
    free(c_part);
    MPI_Type_free(&col_comm);
    MPI_Type_free(&c_type);
    MPI_Finalize();
    return EXIT_SUCCESS;
}