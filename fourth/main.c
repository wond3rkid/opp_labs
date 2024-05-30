#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <float.h>

#define MY_MPI_TAG 111
#define N 256
#define alpha 1e5
#define epsilon 1e-8
int rank, size;
double h_x, h_y, h_z;
//размеры сетки
const int Nx = 512, Ny = 512, Nz = 512;
// границы куба
const double min_x = -1, min_y = -1, min_z = -1;
const double max_x = 1, max_y = 1, max_z = 1;
const double Dx = max_x - min_x, Dy = max_y - min_y, Dz = max_z - min_z;

//  обертка для максимума
double max_abs(double a, double b) {
    return fabs(a) > fabs(b) ? fabs(a) : fabs(b);
}

// вычисление значения фи
double calculate_phi(double x, double y, double z) {
    return x * x + y * y + z * z;
}

// вычисление ро
double calculate_ro(double x, double y, double z) {
    return 6 - alpha * calculate_phi(x, y, z);
}

// обертка для подсчета индекса
int calculate_index(int i, int j, int k) {
    return i * Ny * Nx + j * Nx + k;
}

// проверка, является ли индекс граничным
int is_bound(int index, int N_a) {
    return index == 0 || index == N_a - 1;
}

void initialize_grid(double *grid, int z_part, int start_layer) {
    assert(grid != NULL);
    for (int i = 0; i < z_part + 2; i++) {
        int curr_index = i + start_layer;
        double curr_z = min_z + h_z * curr_index;
        for (int j = 0; j < Nx; j++) {
            double curr_x = min_x + h_x * j;
            for (int k = 0; k < Ny; k++) {
                double curr_y = min_y + h_y * k;
                int index = calculate_index(i, j, k);
                if (is_bound(i, Nz) || is_bound(j, Nx) || is_bound(k, Ny) ||
                    curr_z == min_z || curr_z == max_z) {
                    grid[index] = calculate_phi(curr_x, curr_y, curr_z);
                } else {
                    grid[index] = 0;
                }
            }
        }
    }
}

// обновляем слой линейки и получаем из неё дельту по формула
double update_layer_get_delta(int start_layer, int curr_z, double *grid, double *tmp_grid) {
    int d_index = start_layer + curr_z;
    if (is_bound(d_index, Nz)) {
        memcpy(tmp_grid + curr_z * Nx * Ny, grid + curr_z * Nx * Ny, sizeof(double) * Nx * Ny);
        return 0;
    }
    double delta = 0; // максимальное изменение значения
    double curr_z_coord = min_z + d_index * h_z; // текущая координата по z
    for (int i = 0; i < Nx; i++) {
        double curr_x = min_x + i * h_x;
        for (int j = 0; j < Ny; j++) {
            double curr_y = min_y + j * h_y;
            int index = calculate_index(curr_z, i, j);
            if (is_bound(i, Nx) || is_bound(j, Ny)) {
                tmp_grid[index] = grid[index];
                continue;
            }
            double elem = 0;
            elem += (grid[calculate_index(curr_z, i + 1, j)] + grid[calculate_index(curr_z, i - 1, j)]) / (h_x * h_x);
            elem += (grid[calculate_index(curr_z, i, j + 1)] + grid[calculate_index(curr_z, i, j - 1)]) / (h_y * h_y);
            elem += (grid[calculate_index(curr_z + 1, i, j)] + grid[calculate_index(curr_z - 1, i, j)]) / (h_z * h_z);
            elem -= calculate_ro(curr_x, curr_y, curr_z_coord);
            double factor = 1.0 / (2.0 / (h_x * h_x) + 2.0 / (h_y * h_y) + 2.0 / (h_z * h_z) + alpha);
            tmp_grid[index] = factor * elem;
            delta = max_abs(delta, tmp_grid[index] - grid[index]);
        }
    }
    return delta;
}

int is_solved(double delta) {
    return delta < epsilon ? 1 : 0;
}

// подсчет точности ответа относительно фи
void calculate_accuracy(double *result) {
    double delta = 0;
    for (int i = 0; i < Nz; i++) {
        double z = min_z + i * h_z;
        for (int j = 0; j < Nx; j++) {
            double x = min_x + j * h_x;
            for (int k = 0; k < Ny; k++) {
                double y = min_y + k * h_y;
                double elem = result[calculate_index(i, j, k)];
                double func = calculate_phi(x, y, z);
                delta = max_abs(delta, elem - func);
            }
        }
    }
    if (delta < 100 * epsilon) {
        fprintf(stdout, "Delta is okay\n");
    } else {
        fprintf(stdout, "Delta isn't okay\n");
    }
}

void calculate_jacobi(double *grid, double *tmp_grid, int start_layer, int z_part) {
    double delta = epsilon + 1;
    do {
        double new_delta = 0;
        delta = 0;

        new_delta = update_layer_get_delta(start_layer, 1, grid, tmp_grid);
        delta = max_abs(delta, new_delta);

        new_delta = update_layer_get_delta(start_layer, z_part, grid, tmp_grid);
        delta = max_abs(delta, new_delta);

        MPI_Request requests[4];
        if (rank != 0) {
            MPI_Isend(tmp_grid + Nx * Ny, Nx * Ny, MPI_DOUBLE, rank - 1, MY_MPI_TAG, MPI_COMM_WORLD, &requests[0]);
            MPI_Irecv(tmp_grid, Nx * Ny, MPI_DOUBLE, rank - 1, MY_MPI_TAG, MPI_COMM_WORLD, &requests[1]);
        }
        if (rank != size - 1) {
            MPI_Isend(tmp_grid + z_part * Nx * Ny, Nx * Ny, MPI_DOUBLE, rank + 1, MY_MPI_TAG, MPI_COMM_WORLD,
                      &requests[2]);
            MPI_Irecv(tmp_grid + (z_part + 1) * Nx * Ny, Nx * Ny, MPI_DOUBLE, rank + 1, MY_MPI_TAG, MPI_COMM_WORLD,
                      &requests[3]);
        }

        for (int i = 2; i < z_part; i++) {
            new_delta = update_layer_get_delta(start_layer, i, grid, tmp_grid);
            delta = max_abs(delta, new_delta);
        }

        if (rank != 0) {
            MPI_Wait(&requests[0], MPI_STATUS_IGNORE);
            MPI_Wait(&requests[1], MPI_STATUS_IGNORE);
        }
        if (rank != size - 1) {
            MPI_Wait(&requests[2], MPI_STATUS_IGNORE);
            MPI_Wait(&requests[3], MPI_STATUS_IGNORE);
        }

        memcpy(grid, tmp_grid, (z_part + 2) * Nx * Ny * sizeof(double));
        double max_delta;

        MPI_Allreduce(&delta, &max_delta, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        delta = max_delta;

    } while (!is_solved(delta));

    double *result = NULL;
    if (rank == 0) {
        result = (double *) malloc(sizeof(double) * Nx * Ny * Nz);
    }

    MPI_Gather(grid + Nx * Ny, z_part * Nx * Ny, MPI_DOUBLE, result, z_part * Nx * Ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        calculate_accuracy(result);
        free(result);
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (N % size != 0) {
        if (rank == 0) {
            fprintf(stderr, "Invalid size\n");
        }
        MPI_Finalize();
        exit(1);
    }

    double start = MPI_Wtime();
    h_x = Dx / (Nx - 1);
    h_y = Dy / (Ny - 1);
    h_z = Dz / (Nz - 1);

    int z_part = Nz / size;
    double *grid = (double *) malloc(sizeof(double) * (z_part + 2) * Nx * Ny);
    double *tmp_grid = (double *) malloc(sizeof(double) * (z_part + 2) * Nx * Ny);

    int start_layer = rank * z_part - 1;
    initialize_grid(grid, z_part, start_layer);
    calculate_jacobi(grid, tmp_grid, start_layer, z_part);

    double end = MPI_Wtime();
    if (rank == 0) {
        fprintf(stdout, "Time taken: %lf seconds\n", end - start);
    }

    free(grid);
    free(tmp_grid);
    MPI_Finalize();
    return 0;
}
