#include "sequential_program.h"
#include "parallel_program.h"
#include <time.h>

int main(int argc, char **argv) {
    int N = 1000;
    clock_t start = clock(), end;
    p_preparation_perfomance_free(N);

    end = clock();
    printf("Time taken: %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);
    start = clock();
    preparation_perfomance_free(N);

    end = clock();
    printf("Time taken: %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);
    return 0;
}