#include <mpi.h>
#include <iostream>

#define M 3600
#define N 3600
#define K 3600
#define DIMS 2


void print_matrix(double *matrix, int rows, int columns) {
    if (matrix == nullptr) {
        std::cout << "Unlucky, matrix is null" << std::endl;
    } else {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                std::cout << matrix[i * columns + j] << ' ';
            }
            std::cout << std::endl;
        }
    }
}

void fill_matrix(double *matrix, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix[i * cols + j] = i == j ? 2 : 1;
        }
    }
}

void mult_partly(double *c, const double *a, const double *b, int row_size, int column_size) {
    for (int i = 0; i < row_size; ++i) {
        for (int j = 0; j < column_size; ++j) {
            c[i * column_size + j] = 0;
            for (int k = 0; k < N; ++k) {
                c[i * column_size + j] += a[i * N + k] * b[k * column_size + j];
            }
        }
    }
}

void oneD_comms_create(MPI_Comm comm, MPI_Comm *columns, MPI_Comm *rows) {
    int row_remain_dim[2] = {0, 1};
    int column_remain_dim[2] = {1, 0};
    MPI_Cart_sub(comm, column_remain_dim, columns);
    MPI_Cart_sub(comm, row_remain_dim, rows);
}

void create_mpi_types(MPI_Datatype *b_type, MPI_Datatype *c_type, int row_size, int column_size) {
    // Создание векторов для перемножения без транспонирования
    // количество элементов, размер блока, расстояние между элементами в оригинальном массиве
    MPI_Type_vector(N, column_size, K, MPI_DOUBLE, b_type);
    MPI_Type_vector(row_size, column_size, K, MPI_DOUBLE, c_type);
    MPI_Type_create_resized(*b_type, 0, column_size * sizeof(double), b_type);
    MPI_Type_create_resized(*c_type, 0, column_size * sizeof(double), c_type);
    MPI_Type_commit(b_type);
    MPI_Type_commit(c_type);
}

void fill_vectors(int *b_elem_num, int *b_shifts, const int *dims) {
    //заполнение векторов и сдвигов для дальнейшего умножения
    for (int i = 0; i < dims[1]; ++i) {
        b_shifts[i] = i;
        b_elem_num[i] = 1;
    }
}

void fill_c_vectors(int *c_elem_num, int *c_shifts, int row_size, const int *dims, int grid_size) {
    //сдвиги векторов с
    for (int i = 0; i < grid_size; ++i) {
        c_elem_num[i] = 1;
    }
    for (int i = 0; i < dims[0]; ++i) {
        for (int j = 0; j < dims[1]; ++j) {
            c_shifts[i * dims[1] + j] = i * dims[1] * row_size + j;
        }
    }
}


void calculate(double *A, double *B, double *C, int *dims, int rank, MPI_Comm comm) {
    int row_size = M / dims[0];
    int column_size = K / dims[1];
    int grid_size = 0;
    MPI_Comm_size(comm, &grid_size);
    int coords[2];
    MPI_Cart_coords(comm, rank, DIMS, coords);

    auto *part_a = new double[row_size * N];
    auto *part_b = new double[column_size * N];
    auto *part_c = new double[column_size * row_size];

    int *b_elem_num = nullptr;
    int *b_shifts = nullptr;
    int *c_elem_num = nullptr;
    int *c_shifts = nullptr;
    MPI_Datatype b_type = 0, c_type = 0;
    if (rank == 0) {
        b_elem_num = new int[dims[1]];
        b_shifts = new int[dims[1]];
        c_elem_num = new int[grid_size];
        c_shifts = new int[grid_size];
        create_mpi_types(&b_type, &c_type, row_size, column_size);
        fill_vectors(b_elem_num, b_shifts, dims);
        fill_c_vectors(c_elem_num, c_shifts, row_size, dims, grid_size);
    }

    MPI_Comm column_comm;
    MPI_Comm row_comm;
    oneD_comms_create(comm, &column_comm, &row_comm);
    if (coords[1] == 0) {
        MPI_Scatter(A, row_size * N, MPI_DOUBLE, part_a, row_size * N, MPI_DOUBLE, 0, column_comm);
    }
    if (coords[0] == 0) {
        MPI_Scatterv(B, b_elem_num, b_shifts, b_type, part_b, column_size * N, MPI_DOUBLE, 0, row_comm);
    }
    MPI_Bcast(part_a, row_size * N, MPI_DOUBLE, 0, row_comm);
    MPI_Bcast(part_b, column_size * N, MPI_DOUBLE, 0, column_comm);
    //перемножение кусков матрицы
    mult_partly(part_c, part_a, part_b, row_size, column_size);
    //сбор в нулевом процесе
    MPI_Gatherv(part_c, column_size * row_size, MPI_DOUBLE, C, c_elem_num, c_shifts, c_type, 0, comm);
    if (rank == 0) {
        delete[] b_elem_num;
        delete[] b_shifts;
        delete[] c_elem_num;
        delete[] c_shifts;
        MPI_Type_free(&b_type);
        MPI_Type_free(&c_type);
    }
    delete[] part_a;
    delete[] part_b;
    delete[] part_c;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cout << "Input error" << std::endl;
        return 1;
    }
    // инициализация MPI программы
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // создание 2 мерной сетки размера х на у
    int x, y;
    x = atoi(argv[1]);
    y = atoi(argv[2]);
    int dims[2] = {x, y};
    int periods[2] = {0};
    int reorder = 0;
    MPI_Comm comm;
    MPI_Dims_create(size, DIMS, dims);
    MPI_Cart_create(MPI_COMM_WORLD, DIMS, dims, periods, reorder, &comm);

    double *a = nullptr;
    double *b = nullptr;
    double *c = nullptr;

    if (rank == 0) {
        a = new double[M * N];
        b = new double[K * N];
        c = new double[K * M];
        fill_matrix(a, M, N);
        fill_matrix(b, N, K);
    }
    double startTime = MPI_Wtime();

    calculate(a, b, c, dims, rank, comm);
    if (rank == 0) {
        double endTime = MPI_Wtime();
        //print_matrix(C, M, K);
        std::cout << "Time = " << endTime - startTime << std::endl;
    }
    if (rank == 0) {
        delete[] a;
        delete[] b;
        delete[] c;
    }
    MPI_Finalize();
    return 0;
}
