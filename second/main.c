//#include "sequential_program.h"
#include "parallel_for_program.h"
//#include "parallel_section_program.h"
//#include <time.h>

#define EXIT_ERROR 1
#define EXIT_SUCCESS 0

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Where is N? Mhhhm?\n");
        return EXIT_ERROR;
    }
    int N = atoi(argv[1]);
    fprintf(stdout, "Your N: %d \n", N);
    /* clock_t start = clock(), end;
    preparation_perfomance_free(N);
    end = clock();
    printf("Time taken for non-parallel: %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);
*/
    double start1 = omp_get_wtime();
    p_preparation_perfomance_free(N);
    double end1 = omp_get_wtime();
    printf("Time taken for 'for' parallel: %f seconds\n", end1 - start1);
/*
    double start2 = omp_get_wtime();
    s_preparation_perfomance_free(N);
    double end2 = omp_get_wtime();
    printf("Time taken for 'section' parallel: %f seconds\n", end2 - start2);
  */
    return EXIT_SUCCESS;
}