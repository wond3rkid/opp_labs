#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <mpi.h>
#include <math.h>

#define MAX_TASKS 900
#define ITERATIONS_COUNT 40
#define L 2

int *tasks;
int size;
int rank;
int offset;
double *times;
pthread_mutex_t mutex;

struct job_requester {
    int tasks_count;
    int part_of_res;
};

struct job_requester do_job(int length) {
    int res = 0;
    int tasks_count = 0;
    struct job_requester req;
    while (1) {
        pthread_mutex_lock(&mutex);
        if (offset >= length) {
            pthread_mutex_unlock(&mutex);
            break;
        }
        int current_offset = offset++;
        pthread_mutex_unlock(&mutex);

        int weight = tasks[current_offset];

        for (int j = 0; j < weight; j++) {
            res += (int)sqrt(j);
        }
        tasks_count++;
    }

    req.part_of_res = res;
    req.tasks_count = tasks_count;
    return req;
}

void set_tasks(int iter_count) {
    for (int i = 0; i < size * MAX_TASKS; i++) {
        tasks[i] = abs(50 - i % 100) * abs(rank - iter_count % size) * L;
    }
}

struct job_requester request_tasks() {
    for (int i = 0; i < size; i++) {
        if (i == rank) continue;
        int req_code = 0;
        int help;

        MPI_Send(&req_code, 1, MPI_INT, i, 333, MPI_COMM_WORLD);
        MPI_Recv(&help, 1, MPI_INT, i, 3333, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (help > 0) {
            MPI_Recv(tasks, help, MPI_INT, i, 3333, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            pthread_mutex_lock(&mutex);
            offset = 0;
            pthread_mutex_unlock(&mutex);
        }

        return do_job(help);
    }
}

void find_max_and_min(const double* arr, int len, double* max, double* min) {
    *max = *min = arr[0];
    for (int i = 1; i < len; i++) {
        if (arr[i] > *max) {
            *max = arr[i];
        }
        if (arr[i] < *min) {
            *min = arr[i];
        }
    }
}

void do_tasks() {
    double start = MPI_Wtime();
    struct job_requester res1 = do_job(size * MAX_TASKS);
    struct job_requester res2 = request_tasks();
    double end = MPI_Wtime();

    int tasks_count = res1.tasks_count + res2.tasks_count;
    int part_of_res = res1.part_of_res + res2.part_of_res;
    int res = 0;
    MPI_Allreduce(&part_of_res, &res, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    double work_time = end - start;
    printf("Rank: %d,  tasks done: %d, res = %d, time: %f\n", rank, tasks_count, res, work_time);
    MPI_Gather(&work_time, 1, MPI_DOUBLE, times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double max;
        double min;
        find_max_and_min(times, size, &max, &min);
        double disbalance_time = max - min;
        printf("Disbalance time: %f, proportion: %f\n", disbalance_time, disbalance_time / max * 100);
    }
}

void* work(void* arg) {
    for (int i = 0; i < ITERATIONS_COUNT; i++) {
        pthread_mutex_lock(&mutex);
        offset = 0;
        set_tasks(i);
        pthread_mutex_unlock(&mutex);
        do_tasks();
    }
    int stop_code = 1;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Send(&stop_code, 1, MPI_INT, rank, 333, MPI_COMM_WORLD);
    return NULL;
}

void* listen(void* arg) {
    while (1) {
        int req_code;
        MPI_Status st;

        MPI_Recv(&req_code, 1, MPI_INT, MPI_ANY_SOURCE, 333, MPI_COMM_WORLD, &st);

        if (req_code == 1) {
            break;
        }

        size_t length = size * MAX_TASKS;

        pthread_mutex_lock(&mutex);
        int new_offset = offset + (int)((int)(length - offset) * 0.1);
        int tasks_length = new_offset - offset;
        pthread_mutex_unlock(&mutex);

        MPI_Send(&tasks_length, 1, MPI_INT, st.MPI_SOURCE, 3333, MPI_COMM_WORLD);

        if (tasks_length > 0) {
            pthread_mutex_lock(&mutex);
            int old_offset = offset;
            offset = new_offset;
            pthread_mutex_unlock(&mutex);
            MPI_Send(&tasks[old_offset], tasks_length, MPI_INT, st.MPI_SOURCE, 3333, MPI_COMM_WORLD);
        }
    }
    return NULL;
}

int main(int argc, char* argv[]) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided != MPI_THREAD_MULTIPLE) {
        printf("Invalid thread level support\n");
        exit(-1);
    }

    double start = MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    tasks = (int*) malloc(size * MAX_TASKS * sizeof(int));
    if (rank == 0) {
        times = (double*) malloc(size * sizeof(double));
    }

    pthread_t worker, listener;
    pthread_mutex_init(&mutex, NULL);
    pthread_create(&worker, NULL, work, NULL);
    pthread_create(&listener, NULL, listen, NULL);

    pthread_join(worker, NULL);
    pthread_join(listener, NULL);

    if (tasks != NULL) {
        free(tasks);
    }

    if (rank == 0) {
        free(times);
        double end = MPI_Wtime();
        printf("Time: %f\n", end - start);
    }

    MPI_Finalize();
    return 0;
}
