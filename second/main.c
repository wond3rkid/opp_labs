#include "sequential_program.h"
#include "parallel_program.h"
#include <time.h>

int main(int argc, char **argv) {
    int N = 10000;

    clock_t start = clock(), end;
    preparation_perfomance_free(N);
    end = clock();
    printf("Time taken for non-parallel: %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);

    double start1 = omp_get_wtime();
    p_preparation_perfomance_free(N);
    double end1 = omp_get_wtime();
    printf("Time taken for parallel: %f seconds\n", end1 - start1);
    return 0;
}