#include "sequential_program.h"
#include <time.h>

int main(int argc, char **argv) {
    clock_t start = clock(), end;
    preparation_perfomance_free(1000);
    end = clock();
    printf("Time taken: %f seconds", (double) (end - start) / CLOCKS_PER_SEC);
    return 0;
}