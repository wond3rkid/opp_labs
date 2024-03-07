#include "sequential_program.h"
#include "parallel_program.h"
#include <time.h>

int main(int argc, char **argv) {
    int N = 100;

    clock_t start = clock(), end;
    preparation_perfomance_free(N);
    end = clock();
    printf("Time taken for non-parallel: %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);

    start = clock();
    p_preparation_perfomance_free(N);
    end = clock();
    printf("Time taken for parallel: %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);

    return 0;
}