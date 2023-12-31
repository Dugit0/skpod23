#include "mpi.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define Max(a, b) (((a) > (b)) ? (a) : (b))
#define Min(a, b) (((a) < (b)) ? (a) : (b))
#define  N   10
#define  debug 10

int M;
int m = 5;
double maxeps = 0.1e-7;
int itmax = 100;
/* int i, j, k; */
double check_sum = 0.;
double eps, local_eps;
int rank, num_threads;

double A[N][N][N];


void init() {
    for (int i = 0; i <= N - 1; i++) {
        for (int j = 0; j <= N - 1; j++) {
            for (int k = 0; k <= N - 1; k++) {
                if (i == 0 || i == N - 1 || j == 0 || j == N - 1 || k == 0 || k == N - 1) {
                    A[i][j][k] = 0.;
                } else {
                    A[i][j][k] = (4. + i + j + k);
                }
            }
        }
    }
}

void relax(int st_i, int st_j, int st_k, int tag) {
    int end_i, end_j, end_k;
    end_i = Min(st_i + m, N - 1);
    end_j = Min(st_j + m, N - 1);
    end_k = Min(st_k + m, N - 1);
    st_i = Max(st_i, 1);
    st_j = Max(st_j, 1);
    st_k = Max(st_k, 1);
    if (debug >= 5)
    printf("%d/%d: start = (%d %d %d), end = (%d %d %d), tag = %d\n",
            rank, num_threads, 
            st_i, st_j, st_k, 
            end_i, end_j, end_k, 
            tag);
    /* Вот тут синхронизация с предыдущим */
    if (st_i != 1) {
        int i = st_i - 1;
        int message_len = (end_j - st_j) * (end_k - st_k);
        double *recv_buf = calloc(message_len, sizeof(double));
        if (!recv_buf) {
            printf("Corrupted calloc for recv_buf in %d/%d, message_len = %d\n", rank, num_threads, message_len);
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        }
        int source = ((st_i / m - 1) * (M*M) + (st_j / m) * M + (st_k / m)) % num_threads; // % num_threads???
        /* int tag =  */
        MPI_Status status;
        
        if (debug >= 5)
        printf("%d/%d: reciev from %d, len = %d\n", rank, num_threads, source, message_len);

        MPI_Recv(recv_buf, message_len, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        int send_check = 1;
        MPI_Send(&send_check, 1, MPI_INT, source, tag, MPI_COMM_WORLD);
        for (int j = st_j; j < end_j; j++) {
            for (int k = st_k; k < end_k; k++) {
                A[i][j][k] = recv_buf[(j - st_j)*(end_j - st_j) + (k - st_k)];
            }
        }
        if (debug >= 10)
        printf("Try recv_buf free in %d/%d, addr=%p\n", rank, num_threads, recv_buf);
        free(recv_buf);
        if (debug >= 10)
        printf("Success recv_buf free in %d/%d\n", rank, num_threads);
    }
    if (st_j != 1) {
        int j = st_j - 1;
        int message_len = (end_i - st_i) * (end_k - st_k);
        double *recv_buf = calloc(message_len, sizeof(double));
        if (!recv_buf) {
            printf("Corrupted calloc for recv_buf in %d/%d, message_len = %d\n", rank, num_threads, message_len);
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        }
        int source = ((st_i / m) * (M*M) + (st_j / m - 1) * M + (st_k / m)) % num_threads; // % num_threads???
        /* int tag =  */
        MPI_Status status;
        
        if (debug >= 5)
        printf("%d/%d: reciev from %d, len = %d\n", rank, num_threads, source, message_len);

        MPI_Recv(recv_buf, message_len, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        int send_check = 1;
        MPI_Send(&send_check, 1, MPI_INT, source, tag, MPI_COMM_WORLD);
        for (int i = st_i; i < end_i; i++) {
            for (int k = st_k; k < end_k; k++) {
                A[i][j][k] = recv_buf[(i - st_i)*(end_i - st_i) + (k - st_k)];
            }
        }
        if (debug >= 10)
        printf("Try recv_buf free in %d/%d, addr=%p\n", rank, num_threads, recv_buf);
        free(recv_buf);
        if (debug >= 10)
        printf("Success recv_buf free in %d/%d\n", rank, num_threads);
    }
    if (st_k != 1) {
        int k = st_k - 1;
        int message_len = (end_j - st_j) * (end_i - st_i);
        double *recv_buf = calloc(message_len, sizeof(double));
        if (!recv_buf) {
            printf("Corrupted calloc for recv_buf in %d/%d, message_len = %d\n", rank, num_threads, message_len);
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        }
        int source = ((st_i / m) * (M*M) + (st_j / m) * M + (st_k / m - 1)) % num_threads; // % num_threads???
        /* int tag =  */
        MPI_Status status;
        
        if (debug >= 5)
        printf("%d/%d: reciev from %d, len = %d\n", rank, num_threads, source, message_len);

        MPI_Recv(recv_buf, message_len, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        int send_check = 1;
        MPI_Send(&send_check, 1, MPI_INT, source, tag, MPI_COMM_WORLD);
        for (int i = st_i; i < end_i; i++) {
            for (int j = st_j; j < end_j; j++) {
                A[i][j][k] = recv_buf[(i - st_i)*(end_i - st_i) + (j - st_j)];
            }
        }
        if (debug >= 10)
        printf("Try recv_buf free in %d/%d, addr=%p\n", rank, num_threads, recv_buf);
        free(recv_buf);
        if (debug >= 10)
        printf("Success recv_buf free in %d/%d\n", rank, num_threads);
    }
    for (int i = st_i; i < end_i; i++) { 
        for (int j = st_j; j < end_j; j++) {
            for (int k = st_k; k < end_k; k++) {
                double e;
                e = A[i][j][k];
                A[i][j][k] = (A[i-1][j][k] + A[i+1][j][k] + A[i][j-1][k] + A[i][j+1][k] + A[i][j][k-1] + A[i][j][k+1]) / 6.;
                local_eps = Max(local_eps, fabs(e - A[i][j][k]));
            }    
        }
    }
    /* Вот тут синхронизация с последующим */
    if (end_i != N - 1) {
        int i = end_i - 1;
        int message_len = (end_j - st_j) * (end_k - st_k);
        double *send_buf = calloc(message_len, sizeof(double));
        if (!send_buf) {
            printf("Corrupted calloc for send_buf in %d/%d, message_len = %d\n", rank, num_threads, message_len);
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        }
        for (int j = st_j; j < end_j; j++) {
            for (int k = st_k; k < end_k; k++) {
                send_buf[(j - st_j)*(end_j - st_j) + (k - st_k)] = A[i][j][k];
            }
        }
        int dest = ((end_i / m) * (M*M) + (st_j / m) * M + (st_k / m)) % num_threads; // % num_threads???
        
        if (debug >= 5)
        printf("%d/%d: send to %d, len = %d\n", rank, num_threads, dest, message_len);

        /* MPI_Send(send_buf, message_len, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD); */
        int recv_check = 0;
        MPI_Status status;
        MPI_Sendrecv(send_buf, message_len, MPI_DOUBLE, dest, tag, &recv_check, 1, MPI_INT, dest, tag, MPI_COMM_WORLD, &status);
        if (debug >= 10)
        printf("Try send_buf free in %d/%d, addr=%p\n", rank, num_threads, send_buf);
        free(send_buf);
        if (debug >= 10)
        printf("Success send_buf free in %d/%d\n", rank, num_threads);
    }
    if (end_j != N - 1) {
        int j = end_j - 1;
        int message_len = (end_i - st_i) * (end_k - st_k);
        double *send_buf = calloc(message_len, sizeof(double));
        if (!send_buf) {
            printf("Corrupted calloc for send_buf in %d/%d, message_len = %d\n", rank, num_threads, message_len);
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        }
        for (int i = st_j; i < end_i; i++) {
            for (int k = st_k; k < end_k; k++) {
                send_buf[(i - st_i)*(end_i - st_i) + (k - st_k)] = A[i][j][k];
            }
        }
        int dest = ((st_i / m) * (M*M) + (end_j / m) * M + (st_k / m)) % num_threads; // % num_threads???
        
        if (debug >= 5)
        printf("%d/%d: send to %d, len = %d\n", rank, num_threads, dest, message_len);

        int recv_check = 0;
        MPI_Status status;
        MPI_Sendrecv(send_buf, message_len, MPI_DOUBLE, dest, tag, &recv_check, 1, MPI_INT, dest, tag, MPI_COMM_WORLD, &status);
        if (debug >= 10)
        printf("Try send_buf free in %d/%d, addr=%p\n", rank, num_threads, send_buf);
        free(send_buf);
        if (debug >= 10)
        printf("Success send_buf free in %d/%d\n", rank, num_threads);
    }
    if (end_k != N - 1) {
        int k = end_k - 1;
        int message_len = (end_j - st_j) * (end_i - st_i);
        double *send_buf = calloc(message_len, sizeof(double));
        if (!send_buf) {
            printf("Corrupted calloc for send_buf in %d/%d, message_len = %d\n", rank, num_threads, message_len);
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        }
        
        for (int i = st_i; i < end_i; i++) {
            for (int j = st_j; j < end_j; j++) {
                send_buf[(i - st_i)*(end_i - st_i) + (j - st_j)] = A[i][j][k];
            }
        }
        int dest = ((st_i / m) * (M*M) + (st_j / m) * M + (end_k / m)) % num_threads; // % num_threads???
        
        if (debug >= 5)
        printf("%d/%d: send to %d, len = %d\n", rank, num_threads, dest, message_len);

        int recv_check = 0;
        MPI_Status status;
        MPI_Sendrecv(send_buf, message_len, MPI_DOUBLE, dest, tag, &recv_check, 1, MPI_INT, dest, tag, MPI_COMM_WORLD, &status);
        
        if (debug >= 10)
        printf("Try send_buf free in %d/%d, addr=%p\n", rank, num_threads, send_buf);
        
        free(send_buf);
        
        if (debug >= 10)
        printf("Success send_buf free in %d/%d\n", rank, num_threads);
    }
}


double verify(int st_i, int st_j, int st_k) {
    double s = 0.;
    double local_s = 0.;
    int end_i, end_j, end_k;
    end_i = Min(st_i + m, N);
    end_j = Min(st_j + m, N);
    end_k = Min(st_k + m, N);
    for (int i = st_i; i < end_i; i++) {
        for (int j = st_j; j < end_j; j++) {
            for (int k = st_k; k < end_k; k++) {
                local_s += A[i][j][k] * (i + 1) * (j + 1) * (k + 1) / (N * N * N);
            }
        }
    }
    return local_s;
}


void printa() {
    for (int i = 0; i <= N - 1; i++) {
        for (int j = 0; j <= N - 1; j++) {
            for (int k = 0; k <= N - 1; k++) {
                printf("%f ", A[i][j][k]);
            }
            printf("\n");
        }
        printf("=======================\n");
    }
    return;
}


int main(int argc, char **argv) {

    int status = MPI_Init(&argc, &argv);
    if (status) { 
        printf("MPI not supported\nError with code %d\n", status);
        MPI_Abort(MPI_COMM_WORLD, status);
        return status;
    };

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_threads);
    if (debug >= 5)
    printf("Created %d/%d\n", rank, num_threads);

    M = (N + m - 1) / m;

    double start, end;
    int st_i, st_j, st_k;
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0) {
        start = MPI_Wtime();
    }
    
    init();
    MPI_Barrier(MPI_COMM_WORLD);
    int it;
    /* for (it = 1; it <= itmax; it++) */
    for (it = 1; it <= 1; it++)
    {
        eps = 0.;
        local_eps = 0.;
        for (int cur_cube = rank; cur_cube < M*M*M; cur_cube += num_threads) {
            st_i = cur_cube / (M*M) * m;
            st_j = (cur_cube % (M*M)) / M * m;
            st_k = (cur_cube % M) * m;
            relax(st_i, st_j, st_k, cur_cube / num_threads);
        }
        MPI_Allreduce(&local_eps, &eps, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) {
            printa();
            printf("===============================\n");
            printf("-------------------------------\n");
            printf("===============================\n");
        }
        if (eps < maxeps) {
            break;
        }
    }
    
    double s = 0.;
    for (int cur_cube = rank; cur_cube < M*M*M; cur_cube += num_threads) {
        st_i = cur_cube / (M*M) * m;
        st_j = (cur_cube % (M*M)) / M * m;
        st_k = (cur_cube % M) * m;
        s += verify(st_i, st_j, st_k);
    }
    printf("Proc %d/%d: s = %f\n", rank, num_threads, s);
    MPI_Reduce(&s, &check_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        end = MPI_Wtime();
        printf("check_sum = %lf\n", check_sum);
        printf("Time = %f, num_threads = %d\n", end - start, num_threads);
    }
    MPI_Finalize();

    return 0;
}


