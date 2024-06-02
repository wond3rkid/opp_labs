#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <mpi.h>
#include <math.h>

#define MAX_TASKS 900
#define ITERATIONS_COUNT 40
#define MY_MPI_TAG 123

int *tasks;
int size;
int rank;
int offset;
double *times;
pthread_mutex_t mutex;

struct job_requester {
    int tasks_done;
    int curr_res;
    int weight;
};

// Функция для выполнения заданий
struct job_requester do_job(int length) {
    int res = 0;
    int tasks_done = 0;
    int weight_done = 0;
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
        tasks_done++;
        weight_done += weight;
    }
    req.curr_res = res;
    req.tasks_done = tasks_done;
    req.weight = weight_done;
    return req;
}

// Функция для инициализации задач, установление веса задач
void set_tasks(int iter_count, int rank, int size, int *tasks) {
    int rank_weight = 0;
    for (int i = 0; i < size * MAX_TASKS; i++) {
        tasks[i] = abs(50 - i % 100) * abs(rank - iter_count % size);
        rank_weight += tasks[i];
    }
    fprintf(stdout, "rank : %d | full weight : %d \n", rank, rank_weight);
    int total_weight = 0;
    MPI_Allreduce(&rank_weight, &total_weight, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        fprintf(stdout, "allreduceed total weight : %d \n", total_weight);
    }
}

// Функция для запроса задач у других процессов
struct job_requester request_tasks() {
    struct job_requester final_result = {0, 0, 0};

    for (int i = 0; i < size; i++) {
        if (i == rank) {
            continue;
        }
        int req_code = 0;
        int help;

        MPI_Send(&req_code, 1, MPI_INT, i, MY_MPI_TAG, MPI_COMM_WORLD);
        MPI_Recv(&help, 1, MPI_INT, i, MY_MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (help > 0) {
            MPI_Recv(tasks, help, MPI_INT, i, MY_MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            pthread_mutex_lock(&mutex);
            offset = 0;
            pthread_mutex_unlock(&mutex);

            struct job_requester partial_result = do_job(help);
            final_result.curr_res += partial_result.curr_res;
            final_result.tasks_done += partial_result.tasks_done;
            final_result.weight += partial_result.weight;
            if (partial_result.weight == 0) {
                fprintf(stdout, "proc : %d, weight : %d \n", rank, partial_result.weight);
            }
        }
    }
    if (final_result.weight == 0) {
        fprintf(stdout, "proc : %d, weight : %d \n", rank, final_result.weight);
    }
    return final_result;
}

// Функция для нахождения максимального и минимального значений
void find_max_and_min(const double *array, int len, double *max, double *min) {
    *max = *min = array[0];
    for (int i = 1; i < len; i++) {
        if (array[i] > *max) {
            *max = array[i];
        }
        if (array[i] < *min) {
            *min = array[i];
        }
    }
}

// Основная функция выполнения задач
void do_tasks() {
    double start = MPI_Wtime();
    struct job_requester res1 = do_job(size * MAX_TASKS);
    struct job_requester res2 = request_tasks();
    double end = MPI_Wtime();

    int tasks_done = res1.tasks_done + res2.tasks_done;
    int curr_res = res1.curr_res + res2.curr_res;
    int total_weight = res1.weight + res2.weight;

    int res = 0;
    int all_weights = 0;

    MPI_Allreduce(&curr_res, &res, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&total_weight, &all_weights, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    double work_time = end - start;
    fprintf(stdout, "rank: %d, | tasks done: %d | res = %d | weight: %d | time: %f \n", rank, tasks_done, res, total_weight, work_time);

    MPI_Gather(&work_time, 1, MPI_DOUBLE, times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double max, min;
        find_max_and_min(times, size, &max, &min);
        double disbalance_time = max - min;
        fprintf(stdout, "disbalance time: %f \n", disbalance_time);
        fprintf(stdout, "total weight from all processes: %d\n", all_weights);
    }
}

// Функция потока для выполнения задач
void *work(void *arg) {
    // i 0 0
    for (int i = 0; i < ITERATIONS_COUNT; i++) {
        pthread_mutex_lock(&mutex);
        offset = 0;
        set_tasks(i, rank, size, tasks);
        pthread_mutex_unlock(&mutex);
        do_tasks();
    }
    int stop_code = 1;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Send(&stop_code, 1, MPI_INT, rank, MY_MPI_TAG, MPI_COMM_WORLD);
    return NULL;
}

// Функция потока для прослушивания запросов на задачи
void *listen(void *arg) {
    while (1) {
        int req_code;
        MPI_Status st;

        MPI_Recv(&req_code, 1, MPI_INT, MPI_ANY_SOURCE, MY_MPI_TAG, MPI_COMM_WORLD, &st);

        if (req_code == 1) {
            break;
        }

        size_t length = size * MAX_TASKS;

        pthread_mutex_lock(&mutex);
        int new_offset = offset + (int)((length - offset) * 0.1);
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

int main(int argc, char *argv[]) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided != MPI_THREAD_MULTIPLE) {
        fprintf(stdout, "invalid thread level support\n");
        exit(-1);
    }

    double start = MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    tasks = (int *)malloc(size * MAX_TASKS * sizeof(int));
    if (rank == 0) {
        times = (double *)malloc(size * sizeof(double));
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
        fprintf(stdout, "time: %f\n", end - start);
    }

    MPI_Finalize();
    return 0;
}
