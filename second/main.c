#include "sequential_program.h"
#include "parallel_program.h"
#include <time.h>

int main(int argc, char **argv) {
    clock_t start = clock(), end;
    omp_set_num_threads(omp_get_num_procs());
    preparation_perfomance_free(4);
    p_preparation_perfomance_free(3);

    end = clock();
    printf("Time taken: %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);
    return 0;
}