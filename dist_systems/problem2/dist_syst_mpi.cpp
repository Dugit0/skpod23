#include <mpi.h>
#include <mpi-ext.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstdio>
#include <unistd.h>
#include <signal.h>

#include <string>
#include <limits>

#define Max(a, b) (((a) > (b)) ? (a) : (b))
#define Min(a, b) (((a) < (b)) ? (a) : (b))
#define  debug 0


typedef struct {
    int beg_i;
    int beg_j;
    int beg_k;
    int end_i;
    int end_j;
    int end_k;
} Slice;

typedef struct {
    int i;
    int j;
    int k;
} Triplet;

typedef struct {
    int proc_rank;
    Triplet beg;
    Triplet end;
    Triplet len;
    Triplet send_to;
    Triplet recv_from;
    double ***cube;
} Info;


struct ErrorImitation {
    int flag;
    int rank;
    int iteration;
};

int N;
int M;
int m;
double maxeps = 0.1e-7;
int itmax = 100;
double check_sum = 0.;
double eps = 0;
int rank, num_threads;
Info *metadata;

enum {
    BACKUP_SIZE = 2
};
ErrorImitation breakdowns[BACKUP_SIZE] = {{0, 2, 3}, {0, 100, 7}};
int num_hot_threads;
int break_number = 0;
MPI_Comm comm_world;


Info *init_metadata() {
    Info *res = (Info *) calloc(M*M*M, sizeof(*res));
    for (int i = 0; i < M*M*M; i++) {
        res[i].proc_rank = i % num_hot_threads;
        res[i].beg.i = i / (M*M) * m;
        res[i].beg.j = i % (M*M) / M * m;
        res[i].beg.k = i % M * m;
        res[i].end.i = Min(res[i].beg.i + m - 1, N - 1);
        res[i].end.j = Min(res[i].beg.j + m - 1, N - 1);
        res[i].end.k = Min(res[i].beg.k + m - 1, N - 1);
        res[i].len.i = res[i].end.i - res[i].beg.i + 1;
        res[i].len.j = res[i].end.j - res[i].beg.j + 1;
        res[i].len.k = res[i].end.k - res[i].beg.k + 1;

        int buf;
        // Send to
        buf = (i / (M*M) + 1) * (M*M) + (i % (M*M));
        if (i / (M*M) + 1 >= M) {
            buf = -1;
        }
        res[i].send_to.i = buf;

        buf = i / (M*M) * (M*M) + (i % (M*M) / M + 1) * M + (i % M);
        if (i % (M*M) / M + 1 >= M) {
            buf = -1;
        }
        res[i].send_to.j = buf;

        buf = (i / M * M) + (i % M + 1);
        if (i % M + 1 >= M) {
            buf = -1;
        }
        res[i].send_to.k = buf;

        // Recv from
        buf = (i / (M*M) - 1) * (M*M) + (i % (M*M));
        if (i / (M*M) - 1 < 0) {
            buf = -1;
        }
        res[i].recv_from.i = buf;

        buf = i / (M*M) * (M*M) + (i % (M*M) / M - 1) * M + (i % M);
        if (i % (M*M) / M - 1 < 0) {
            buf = -1;
        }
        res[i].recv_from.j = buf;

        buf = (i / M * M) + (i % M - 1);
        if (i % M - 1 < 0) {
            buf = -1;
        }
        res[i].recv_from.k = buf;
    }

    return res;
}


double ***init(Info info) {

    printf("%d/%d: Start init\n", rank, num_threads);
    double ***res = (double ***) calloc(info.len.i, sizeof(*res));

    for (int i = 0; i < info.len.i; i++) {
        res[i] = (double **) calloc(info.len.j, sizeof(*(res[i])));
        for (int j = 0; j < info.len.j; j++) {
            res[i][j] = (double *) calloc(info.len.k, sizeof(*(res[i][j])));
            for (int k = 0; k < info.len.k; k++) {
                int real_i = info.beg.i + i;
                int real_j = info.beg.j + j;
                int real_k = info.beg.k + k;
                if (real_i == 0 || real_i == N - 1 || real_j == 0 || real_j == N - 1 || real_k == 0 || real_k == N - 1) {
                    res[i][j][k] = 0.;
                } else {
                    res[i][j][k] = (4. + real_i + real_j + real_k);
                }
                // printf("%f\n", res[i][j][k]);
            }
        }
    }
    printf("%d/%d: Finish init\n", rank, num_threads);
    return res;
}


void free_cube(Info info) {
    for (int i = 0; i < info.len.i; i++) {
        for (int j = 0; j < info.len.j; j++) {
            free(info.cube[i][j]);
        }
        free(info.cube[i]);
    }
    free(info.cube);
}


double relax(Info info) {

    // Отправка в предыдущий для вычислений в нем
    // по оси i
    if (info.recv_from.i != -1) {
        if (debug >= 5)
            printf("%d/%d: Send to prev for calc i\n", rank, num_threads);
        int message_len = info.len.j * info.len.k;
        double *send_buf = (double *) calloc(message_len, sizeof(*send_buf));
        int i = 0;
        for (int j = 0; j < info.len.j; j++) {
            for (int k = 0; k < info.len.k; k++) {
                int ind = j * info.len.k + k;
                send_buf[ind] = info.cube[i][j][k];
            }
        }
        int recv_check = 0;
        MPI_Status status;
        MPI_Sendrecv(send_buf, message_len, MPI_DOUBLE, info.recv_from.i, 0,
                &recv_check, 1, MPI_INT, info.recv_from.i, 0,
                comm_world, &status);
        free(send_buf);
    }
    // по оси j
    if (info.recv_from.j != -1) {
        if (debug >= 5)
            printf("%d/%d: Send to prev for calc j\n", rank, num_threads);
        int message_len = info.len.i * info.len.k;
        double *send_buf = (double *) calloc(message_len, sizeof(*send_buf));
        int j = 0;
        for (int i = 0; i < info.len.i; i++) {
            for (int k = 0; k < info.len.k; k++) {
                int ind = i * info.len.k + k;
                send_buf[ind] = info.cube[i][j][k];
            }
        }
        int recv_check = 0;
        MPI_Status status;
        MPI_Sendrecv(send_buf, message_len, MPI_DOUBLE, info.recv_from.j, 0,
                &recv_check, 1, MPI_INT, info.recv_from.j, 0,
                comm_world, &status);
        free(send_buf);
    }
    // по оси k
    if (info.recv_from.k != -1) {
        if (debug >= 5)
            printf("%d/%d: Send to prev for calc k\n", rank, num_threads);
        int message_len = info.len.i * info.len.j;
        double *send_buf = (double *) calloc(message_len, sizeof(*send_buf));
        int k = 0;
        for (int i = 0; i < info.len.i; i++) {
            for (int j = 0; j < info.len.j; j++) {
                int ind = i * info.len.j + j;
                send_buf[ind] = info.cube[i][j][k];
            }
        }
        int recv_check = 0;
        MPI_Status status;
        MPI_Sendrecv(send_buf, message_len, MPI_DOUBLE, info.recv_from.k, 0,
                &recv_check, 1, MPI_INT, info.recv_from.k, 0,
                comm_world, &status);
        free(send_buf);
    }

    // Формирование рабочего буфера
    int real_len_i = (info.recv_from.i == -1) ? info.len.i : (info.len.i + 1);
    int real_len_j = (info.recv_from.j == -1) ? info.len.j : (info.len.j + 1);
    int real_len_k = (info.recv_from.k == -1) ? info.len.k : (info.len.k + 1);

    real_len_i = (info.send_to.i == -1) ? real_len_i : (real_len_i + 1);
    real_len_j = (info.send_to.j == -1) ? real_len_j : (real_len_j + 1);
    real_len_k = (info.send_to.k == -1) ? real_len_k : (real_len_k + 1);

    // Выделение рабочего буфера
    if (debug >= 5)
        printf("%d/%d: Alloc work buf\n", rank, num_threads);

    double ***work_buf = (double ***) calloc(real_len_i, sizeof(*work_buf));
    for (int i = 0; i < real_len_i; i++) {
        work_buf[i] = (double **) calloc(real_len_j, sizeof(*(work_buf[i])));
        for (int j = 0; j < real_len_j; j++) {
            work_buf[i][j] = (double *) calloc(real_len_k, sizeof(*(work_buf[i][j])));
        }
    }

    // Получение из последующих невычисленных
    // по оси k
    if (debug >= 5)
        printf("%d/%d: Recv from next, what not calc k\n", rank, num_threads);
    {
        int message_len = info.len.i * info.len.j;
        double *recv_buf = (double *) calloc(message_len, sizeof(*recv_buf));
        if (info.send_to.k == -1) {
            for (int ind = 0; ind < message_len; ind++) {
                recv_buf[ind] = 0;
            }
        } else {
            MPI_Status status;
            MPI_Recv(recv_buf, message_len, MPI_DOUBLE,
                    info.send_to.k, 0, comm_world, &status);
            int send_check = 1;
            MPI_Send(&send_check, 1, MPI_INT, info.send_to.k, 0, comm_world);
        }
        int k = real_len_k - 1;
        for (int i = 0; i < info.len.i; i++) {
            for (int j = 0; j < info.len.j; j++) {
                int ind = i * info.len.j + j;
                work_buf[i][j][k] = recv_buf[ind];
            }
        }
        free(recv_buf);
    }
    // по оси j
    if (debug >= 5)
        printf("%d/%d: Recv from next, what not calc j\n", rank, num_threads);
    {
        int message_len = info.len.i * info.len.k;
        double *recv_buf = (double *) calloc(message_len, sizeof(*recv_buf));
        if (info.send_to.j == -1) {
            for (int ind = 0; ind < message_len; ind++) {
                recv_buf[ind] = 0;
            }
        } else {
            MPI_Status status;
            MPI_Recv(recv_buf, message_len, MPI_DOUBLE,
                    info.send_to.j, 0, comm_world, &status);
            int send_check = 1;
            MPI_Send(&send_check, 1, MPI_INT, info.send_to.j, 0, comm_world);
        }
        int j = real_len_j - 1;
        for (int i = 0; i < info.len.i; i++) {
            for (int k = 0; k < info.len.k; k++) {
                int ind = i * info.len.k + k;
                work_buf[i][j][k] = recv_buf[ind];
            }
        }
        free(recv_buf);
    }
    // по оси i
    if (debug >= 5)
        printf("%d/%d: Recv from next, what not calc i\n", rank, num_threads);
    {
        int message_len = info.len.j * info.len.k;
        double *recv_buf = (double *) calloc(message_len, sizeof(*recv_buf));
        if (info.send_to.i == -1) {
            for (int ind = 0; ind < message_len; ind++) {
                recv_buf[ind] = 0;
            }
        } else {
            MPI_Status status;
            MPI_Recv(recv_buf, message_len, MPI_DOUBLE,
                    info.send_to.i, 0, comm_world, &status);
            int send_check = 1;
            MPI_Send(&send_check, 1, MPI_INT, info.send_to.i, 0, comm_world);
        }
        int i = real_len_i - 1;
        for (int j = 0; j < info.len.j; j++) {
            for (int k = 0; k < info.len.k; k++) {
                int ind = j * info.len.k + k;
                work_buf[i][j][k] = recv_buf[ind];
            }
        }
        free(recv_buf);
    }


    // Получение из предыдущих вычисленных
    // по оси i
    if (debug >= 5)
        printf("%d/%d: Recv from prev, what calc i\n", rank, num_threads);
    {
        int message_len = info.len.j * info.len.k;
        double *recv_buf = (double *) calloc(message_len, sizeof(*recv_buf));
        if (info.recv_from.i == -1) {
            for (int ind = 0; ind < message_len; ind++) {
                recv_buf[ind] = 0;
            }
        } else {
            MPI_Status status;
            MPI_Recv(recv_buf, message_len, MPI_DOUBLE,
                    info.recv_from.i, 0, comm_world, &status);
            int send_check = 1;
            MPI_Send(&send_check, 1, MPI_INT, info.recv_from.i, 0, comm_world);
        }
        int i = 0;
        for (int j = 0; j < info.len.j; j++) {
            for (int k = 0; k < info.len.k; k++) {
                int ind = j * info.len.k + k;
                work_buf[i][j][k] = recv_buf[ind];
            }
        }
        free(recv_buf);
    }
    // по оси j
    if (debug >= 5)
        printf("%d/%d: Recv from prev, what calc j\n", rank, num_threads);
    {
        int message_len = info.len.i * info.len.k;
        double *recv_buf = (double *) calloc(message_len, sizeof(*recv_buf));
        if (info.recv_from.j == -1) {
            for (int ind = 0; ind < message_len; ind++) {
                recv_buf[ind] = 0;
            }
        } else {
            MPI_Status status;
            MPI_Recv(recv_buf, message_len, MPI_DOUBLE,
                    info.recv_from.j, 0, comm_world, &status);
            int send_check = 1;
            MPI_Send(&send_check, 1, MPI_INT, info.recv_from.j, 0, comm_world);
        }
        int j = 0;
        for (int i = 0; i < info.len.i; i++) {
            for (int k = 0; k < info.len.k; k++) {
                int ind = i * info.len.k + k;
                work_buf[i][j][k] = recv_buf[ind];
            }
        }
        free(recv_buf);
    }
    // по оси k
    if (debug >= 5)
        printf("%d/%d: Recv from prev, what calc k\n", rank, num_threads);
    {
        int message_len = info.len.i * info.len.j;
        double *recv_buf = (double *) calloc(message_len, sizeof(*recv_buf));
        if (info.recv_from.k == -1) {
            for (int ind = 0; ind < message_len; ind++) {
                recv_buf[ind] = 0;
            }
        } else {
            MPI_Status status;
            MPI_Recv(recv_buf, message_len, MPI_DOUBLE,
                    info.recv_from.k, 0, comm_world, &status);
            int send_check = 1;
            MPI_Send(&send_check, 1, MPI_INT, info.recv_from.k, 0, comm_world);
        }
        int k = 0;
        for (int i = 0; i < info.len.i; i++) {
            for (int j = 0; j < info.len.j; j++) {
                int ind = i * info.len.j + j;
                work_buf[i][j][k] = recv_buf[ind];
            }
        }
        free(recv_buf);
    }


    // Копирование в рабочий буфер
    if (debug >= 5)
        printf("%d/%d: Copy to work buf\n", rank, num_threads);
    for (int i = 0; i < info.len.i; i++) {
        for (int j = 0; j < info.len.j; j++) {
            for (int k = 0; k < info.len.k; k++) {
                int shift_i = (info.recv_from.i == -1) ? 0 : 1;
                int shift_j = (info.recv_from.j == -1) ? 0 : 1;
                int shift_k = (info.recv_from.k == -1) ? 0 : 1;
                int work_i = i + shift_i;
                int work_j = j + shift_j;
                int work_k = k + shift_k;
                work_buf[work_i][work_j][work_k] = info.cube[i][j][k];
            }
        }
    }

    // Релаксация
    if (debug >= 5)
        printf("%d/%d: Relax\n", rank, num_threads);
    double local_eps = 0.;
    for (int i = 1; i < real_len_i - 1; i++) {
        for (int j = 1; j < real_len_j - 1; j++) {
            for (int k = 1; k < real_len_k - 1; k++) {
                double e;
                e = work_buf[i][j][k];
                work_buf[i][j][k] = (work_buf[i-1][j][k] + work_buf[i+1][j][k] + work_buf[i][j-1][k] + work_buf[i][j+1][k] + work_buf[i][j][k-1] + work_buf[i][j][k+1]) / 6.;
                local_eps = Max(local_eps, fabs(e - work_buf[i][j][k]));
            }
        }
    }

    // Копирование из рабочего буфера
    if (debug >= 5)
        printf("%d/%d: Copy from work buf\n", rank, num_threads);
    for (int i = 0; i < info.len.i; i++) {
        for (int j = 0; j < info.len.j; j++) {
            for (int k = 0; k < info.len.k; k++) {
                int shift_i = (info.recv_from.i == -1) ? 0 : 1;
                int shift_j = (info.recv_from.j == -1) ? 0 : 1;
                int shift_k = (info.recv_from.k == -1) ? 0 : 1;
                int work_i = i + shift_i;
                int work_j = j + shift_j;
                int work_k = k + shift_k;
                info.cube[i][j][k] = work_buf[work_i][work_j][work_k];
            }
        }
    }

    // Освобождение рабочего буфера
    if (debug >= 5)
        printf("%d/%d: Free work_buf\n", rank, num_threads);
    for (int i = 0; i < real_len_i; i++) {
        for (int j = 0; j < real_len_j; j++) {
            free(work_buf[i][j]);
        }
        free(work_buf[i]);
    }
    free(work_buf);

    // Отправка в последующий вычисленного
    // по оси k
    if (info.send_to.k != -1) {
        if (debug >= 5)
            printf("%d/%d: Send to next k\n", rank, num_threads);
        int message_len = info.len.i * info.len.j;
        double *send_buf = (double *) calloc(message_len, sizeof(*send_buf));
        int k = info.end.k;
        for (int i = 0; i < info.len.i; i++) {
            for (int j = 0; j < info.len.j; j++) {
                int ind = i * info.len.j + j;
                send_buf[ind] = info.cube[i][j][k];
            }
        }
        int recv_check = 0;
        MPI_Status status;
        MPI_Sendrecv(send_buf, message_len, MPI_DOUBLE, info.send_to.k, 0,
                &recv_check, 1, MPI_INT, info.send_to.k, 0,
                comm_world, &status);
        free(send_buf);
    }
    // по оси j
    if (info.send_to.j != -1) {
        if (debug >= 5)
            printf("%d/%d: Send to next j\n", rank, num_threads);
        int message_len = info.len.i * info.len.k;
        double *send_buf = (double *) calloc(message_len, sizeof(*send_buf));
        int j = info.end.j;
        for (int i = 0; i < info.len.i; i++) {
            for (int k = 0; k < info.len.k; k++) {
                int ind = i * info.len.k + k;
                send_buf[ind] = info.cube[i][j][k];
            }
        }
        int recv_check = 0;
        MPI_Status status;
        MPI_Sendrecv(send_buf, message_len, MPI_DOUBLE, info.send_to.j, 0,
                &recv_check, 1, MPI_INT, info.send_to.j, 0,
                comm_world, &status);
        free(send_buf);
    }
    // по оси i
    if (info.send_to.i != -1) {
        if (debug >= 5)
            printf("%d/%d: Send to next i\n", rank, num_threads);
        int message_len = info.len.j * info.len.k;
        double *send_buf = (double *) calloc(message_len, sizeof(*send_buf));
        int i = info.end.i;
        for (int j = 0; j < info.len.j; j++) {
            for (int k = 0; k < info.len.k; k++) {
                int ind = j * info.len.k + k;
                send_buf[ind] = info.cube[i][j][k];
            }
        }
        int recv_check = 0;
        MPI_Status status;
        MPI_Sendrecv(send_buf, message_len, MPI_DOUBLE, info.send_to.i, 0,
                &recv_check, 1, MPI_INT, info.send_to.i, 0,
                comm_world, &status);
        free(send_buf);
    }

    return local_eps;
}

// void relax()
// {

//     for (i = 1; i <= N - 2; i++) {
//         for (j = 1; j <= N - 2; j++) {
//             for (k = 1; k <= N - 2; k++) {
//                 double e;
//                 e = A[i][j][k];
//                 A[i][j][k] = (A[i-1][j][k] + A[i+1][j][k] + A[i][j-1][k] + A[i][j+1][k] + A[i][j][k-1] + A[i][j][k+1]) / 6.;
//                 eps=Max(eps, fabs(e - A[i][j][k]));
//             }
//         }
//     }
// }



double verify(Info info) {
    double local_s = 0.;
    for (int i = 0; i < info.len.i; i++) {
        for (int j = 0; j < info.len.j; j++) {
            for (int k = 0; k < info.len.k; k++) {
                int real_i = info.beg.i + i;
                int real_j = info.beg.j + j;
                int real_k = info.beg.k + k;
                local_s += info.cube[i][j][k] * (real_i + 1) * (real_j + 1) * (real_k + 1) / (N * N * N);
                // printf("%f\n", local_s);
            }
        }
    }
    return local_s;
}


void print_cube(Info info);
void print_info(Info i);
void print_metadata(Info *metadata);


static void error_handler_function(MPI_Comm*, int*, ...) {
    breakdowns[break_number].flag = 1;
    break_number++;
    int old_rank = rank;

    MPIX_Comm_revoke(comm_world);
    MPIX_Comm_shrink(comm_world, &comm_world);

    MPI_Comm_rank(comm_world, &rank);
    MPI_Comm_size(comm_world, &num_threads);
    printf("%d/%d: Breakdown %d. old_rank = %d, new_rank = %d\n",
           rank, num_threads, break_number,old_rank, rank);
    usleep(100);

    throw 0;
}


void make_backup(Info info) {
    std::string filename = "backups/" + std::to_string(rank) + ".mpi";
    std::FILE* file = std::fopen(filename.c_str(), "w");
    for (int i = 0; i < info.len.i; i++) {
        for (int j = 0; j < info.len.j; j++) {
            for (int k = 0; k < info.len.k; k++) {
                fprintf(file, "%lf ", info.cube[i][j][k]);
            }
        }
    }
    std::fclose(file);
    // MPI_File file;
    // MPI_File_open(MPI_COMM_SELF, filename.c_str(),
    //         MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    // for (int cur_cube = rank; cur_cube < M*M*M; cur_cube += num_hot_threads) {
    //     Info& info = metadata[cur_cube];
    //     for (int i = 0; i < info.len.i; i++) {
    //         for (int j = 0; j < info.len.j; j++) {
    //             for (int k = 0; k < info.len.k; k++) {
    //                 MPI_Status status;
    //                 MPI_File_write(file, &info.cube[i][j][k], 1, MPI_DOUBLE, &status);
    //             }
    //         }
    //     }
    // }
    // MPI_File_close(&file);
    return;
}


// void load_backup(Info* metadata) {


void load_backup() {
    printf("%d/%d: Try load backup\n", rank, num_threads);
    printf("%d/%d: metadata = %p\n", rank, num_threads, metadata);


    // std::string filename = "backups/" + std::to_string(rank) + ".mpi";
    // std::FILE* file = std::fopen(filename.c_str(), "r");
    // // printf("%d/%d: Open backup file %p\n", rank, num_threads, file);
    // for (int cur_cube = rank; cur_cube < M*M*M; cur_cube += num_hot_threads) {
    //     printf("%d/%d: Try get metadata = %p\n", rank, num_threads, metadata);
    //     Info info = metadata[cur_cube];
    //     printf("%d/%d: Try init cube\n", rank, num_threads);
    //     info.cube = init(info);
    //     printf("%d/%d: Init cube\n", rank, num_threads);
    //     for (int i = 0; i < info.len.i; i++) {
    //         for (int j = 0; j < info.len.j; j++) {
    //             for (int k = 0; k < info.len.k; k++) {
    //                 // printf("%d/%d: Try write %d %d %d\n",
    //                 //         rank, num_threads, i, j, k);
    //                 info.cube[i][j][k] = 0;
    //                 // fscanf(file, "%lf ", &(info.cube[i][j][k]));
    //             }
    //         }
    //     }
    //     printf("%d/%d: Load cube\n", rank, num_threads);
    // }
    // std::fclose(file);


    // MPI_File file;
    // MPI_File_open(MPI_COMM_SELF, filename.c_str(),
    //         MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    // for (int cur_cube = rank; cur_cube < M*M*M; cur_cube += num_hot_threads) {
    //     Info& info = metadata[cur_cube];
    //     for (int i = 0; i < info.len.i; i++) {
    //         for (int j = 0; j < info.len.j; j++) {
    //             for (int k = 0; k < info.len.k; k++) {
    //                 MPI_Status status;
    //                 MPI_File_read(file, &info.cube[i][j][k], 1, MPI_DOUBLE, &status);
    //             }
    //         }
    //     }
    // }
    // MPI_File_close(&file);
    return;
}



int main(int argc, char **argv) {

    int status = MPI_Init(&argc, &argv);
    if (status) {
        printf("MPI not supported\nError with code %d\n", status);
        MPI_Abort(comm_world, status);
        return status;
    };
    comm_world = MPI_COMM_WORLD;
    MPI_Comm_rank(comm_world, &rank);
    MPI_Comm_size(comm_world, &num_threads);
    num_hot_threads = num_threads - BACKUP_SIZE;

    // Error handler
    MPI_Errhandler error_handler;
    MPI_Comm_create_errhandler(error_handler_function, &error_handler);
    MPI_Comm_set_errhandler(comm_world, error_handler);

    if (debug >= 5)
        printf("Created %d/%d\n", rank, num_threads);
    MPI_Barrier(comm_world);

    N = 120;
    // Нужно вычислять m(N, num_hot_threads)!!!
    m = 60;
    M = (N + m - 1) / m;


    if (rank == 0) {
        if (debug >= 1)
            printf("N = %d, m = %d, M = %d\n", N, m, M);
    }
    Info *metadata = init_metadata();

    if (rank == 0) {
        if (debug >= 1)
            print_metadata(metadata);
    }
    MPI_Barrier(comm_world);

    double start, end;
    if (rank == 0) {
        start = MPI_Wtime();
    }

    int len_cubes = (M*M*M - rank) / num_hot_threads + 1;
    int *nums_cubes = (int *) calloc(len_cubes, sizeof(*nums_cubes));  // Зачем???

    int tmp_ind;
    tmp_ind = 0;
    for (int cur_cube = rank; cur_cube < M*M*M; cur_cube += num_hot_threads) {
        nums_cubes[tmp_ind] = cur_cube;
        metadata[cur_cube].cube = init(metadata[cur_cube]);
    }

    MPI_Barrier(comm_world);

    for (int it = 1; it <= itmax; it++) {
        usleep(1000);
        try {
            if (!breakdowns[break_number].flag
                    && rank == breakdowns[break_number].rank
                    && it == breakdowns[break_number].iteration) {
                raise(SIGKILL);
            }
            if (rank < num_hot_threads) {
                eps = 0.;
                double proc_eps = 0.;
                for (int cur_cube = rank; cur_cube < M*M*M; cur_cube += num_hot_threads) {
                    double cur_cube_eps = relax(metadata[cur_cube]);
                    proc_eps = Max(proc_eps, cur_cube_eps);
                    make_backup(metadata[cur_cube]);
                }
                MPI_Barrier(comm_world);
                MPI_Allreduce(&proc_eps, &eps, 1, MPI_DOUBLE, MPI_MAX, comm_world);
                MPI_Barrier(comm_world);
                // if (rank == 0) {
                //     printf("eps = %f\n", eps);
                // }
                if (eps < maxeps) {
                    break;
                }
            } else {
                MPI_Barrier(comm_world);
                eps = 0.;
                double proc_eps = std::numeric_limits<double>::min();
                MPI_Allreduce(&proc_eps, &eps, 1, MPI_DOUBLE, MPI_MAX, comm_world);
                MPI_Barrier(comm_world);
            }
        } catch (int error) {
            it--;
            printf("%d/%d: metadata = %p\n", rank, num_threads, metadata);
            load_backup();
            MPI_Barrier(comm_world);
        }
    }

    double s = 0.;
    for (int cur_cube = rank; cur_cube < M*M*M; cur_cube += num_hot_threads) {
        s += verify(metadata[cur_cube]);
    }
    // printf("Proc %d/%d: s = %f\n", rank, num_threads, s);
    MPI_Reduce(&s, &check_sum, 1, MPI_DOUBLE, MPI_SUM, 0, comm_world);
    MPI_Barrier(comm_world);
    if (rank == 0) {
        end = MPI_Wtime();
        printf("eps = %f\n", eps);
        printf("check_sum = %lf\n", check_sum);
        printf("Time = %f, num_threads = %d\n", end - start, num_threads);
    }
    for (int cur_cube = rank; cur_cube < M*M*M; cur_cube += num_hot_threads) {
        free_cube(metadata[cur_cube]);
    }
    free(nums_cubes);

    free(metadata);
    MPI_Finalize();

    return 0;
}



// Debug print functions

void print_cube(Info info) {
    print_info(info);
    for (int i = 0; i < info.len.i; i++) {
        for (int j = 0; j < info.len.j; j++) {
            for (int k = 0; k < info.len.k; k++) {
                printf("%f ", info.cube[i][j][k]);
            }
            printf("\n");
        }
        printf("=======================\n");
    }
}


void print_info(Info i) {
    printf("rank = %d, begin = (%d %d %d), end = (%d %d %d)\nsend_to = (%d %d %d), recv_from = (%d %d %d)\n",
            i.proc_rank,
            i.beg.i, i.beg.j, i.beg.k,
            i.end.i, i.end.j, i.end.k,
            i.send_to.i, i.send_to.j, i.send_to.k,
            i.recv_from.i, i.recv_from.j, i.recv_from.k);
}


void print_metadata(Info *metadata) {

    printf("Struct:\n");
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < M; k++) {
                int ind = i * M * M + j * M + k;
                printf("%2d | ", ind);
            }
            printf("\n");
            // printf("\n--------------\n");
        }
        printf("===========\n");
    }

    printf("\nCubes:\n");
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < M; k++) {
                int ind = i * M * M + j * M + k;
                printf("%2d %2d %2d | ", metadata[ind].beg.i, metadata[ind].beg.j, metadata[ind].beg.k);
            }
            printf("\n");
            for (int k = 0; k < M; k++) {
                int ind = i * M * M + j * M + k;
                printf("%2d %2d %2d | ", metadata[ind].end.i, metadata[ind].end.j, metadata[ind].end.k);
            }
            printf("\n");
            for (int k = 0; k < M; k++) {
                int ind = i * M * M + j * M + k;
                printf("%2d %2d %2d | ", metadata[ind].len.i, metadata[ind].len.j, metadata[ind].len.k);
            }
            printf("\n");
            printf("-------------\n");
        }
        printf("===========\n");
    }

    printf("\nSend to:\n");
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < M; k++) {
                int ind = i * M * M + j * M + k;
                printf("%2d %2d %2d | ", metadata[ind].send_to.i, metadata[ind].send_to.j, metadata[ind].send_to.k);
            }
            printf("\n");
        }
        printf("===========\n");
    }

    printf("\nRecv from:\n");
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < M; k++) {
                int ind = i * M * M + j * M + k;
                printf("%2d %2d %2d | ", metadata[ind].recv_from.i, metadata[ind].recv_from.j, metadata[ind].recv_from.k);
            }
            printf("\n");
        }
        printf("===========\n");
    }
}
