#include "sequential_program.h"
#include "parallel_program.h"
#include <time.h>

int main(int argc, char **argv) {
    clock_t start = clock(), end;
    printf("Initial number of threads: %d\n", omp_get_num_threads());
    omp_set_num_threads(4);
    printf("Number of threads after setting: %d\n", omp_get_num_threads());
    p_preparation_perfomance_free(10000);
    end = clock();
    printf("Time taken: %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);
    return 0;
}