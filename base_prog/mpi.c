#include "mpi.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define Max(a, b) (((a) > (b)) ? (a) : (b))
#define Min(a, b) (((a) < (b)) ? (a) : (b))
#define  N   10

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
    /* Вот тут синхронизация с предыдущим */
    if (st_i != 1) {
        int i = st_i - 1;
        int message_len = (end_j - st_j) * (end_k - st_k);
        double *recv_buf = calloc(message_len, sizeof(double));
        int source = ((st_i / m - 1) * (M*M) + (st_j / m) * M + (st_k / m)) % num_threads; // % num_threads???
        /* int tag =  */
        MPI_Status status;
        MPI_Recv(recv_buf, message_len, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        for (int j = st_j; j < end_j; j++) {
            for (int k = st_k; k < end_k; k++) {
                A[i][j][k] = recv_buf[j*N + k];
            }
        }
    }
    if (st_j != 1) {
        int j = st_j - 1;
        int message_len = (end_i - st_i) * (end_k - st_k);
        double *recv_buf = calloc(message_len, sizeof(double));
        int source = ((st_i / m) * (M*M) + (st_j / m - 1) * M + (st_k / m)) % num_threads; // % num_threads???
        /* int tag =  */
        MPI_Status status;
        MPI_Recv(recv_buf, message_len, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        for (int i = st_i; i < end_i; i++) {
            for (int k = st_k; k < end_k; k++) {
                A[i][j][k] = recv_buf[i*N + k];
            }
        }
    }
    if (st_k != 1) {
        int k = st_k - 1;
        int message_len = (end_j - st_j) * (end_i - st_i);
        double *recv_buf = calloc(message_len, sizeof(double));
        int source = ((st_i / m) * (M*M) + (st_j / m) * M + (st_k / m - 1)) % num_threads; // % num_threads???
        /* int tag =  */
        MPI_Status status;
        MPI_Recv(recv_buf, message_len, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        for (int i = st_i; i < end_i; i++) {
            for (int j = st_j; j < end_j; j++) {
                A[i][j][k] = recv_buf[i*N + j];
            }
        }
    }
    for (int i = st_i; i < end_i; i++) { 
        for (int j = st_j; j < end_j; j++) {
            for (int k = st_k; k < end_k; k++) {
                double e;
                e = A[i][j][k];
                A[i][j][k] = (A[i-1][j][k] + A[i+1][j][k] + A[i][j-1][k] + A[i][j+1][k] + A[i][j][k-1] + A[i][j][k+1]) / 6.;
                local_eps=Max(local_eps, fabs(e - A[i][j][k]));
            }    
        }
    }
    /* Вот тут синхронизация с последующим */
    if (end_i != N - 1) {
        int i = end_i - 1;
        int message_len = (end_j - st_j) * (end_k - st_k);
        double *send_buf = calloc(message_len, sizeof(double));
        for (int j = st_j; j < end_j; j++) {
            for (int k = st_k; k < end_k; k++) {
                send_buf[j*N + k] = A[i][j][k];
            }
        }
        int dest = ((end_i / m) * (M*M) + (st_j / m) * M + (st_k / m)) % num_threads; // % num_threads???
        MPI_Send(send_buf, message_len, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
    }
    if (end_j != N - 1) {
        int j = end_j - 1;
        int message_len = (end_i - st_i) * (end_k - st_k);
        double *send_buf = calloc(message_len, sizeof(double));
        for (int i = st_j; i < end_i; i++) {
            for (int k = st_k; k < end_k; k++) {
                send_buf[i*N + k] = A[i][j][k];
            }
        }
        int dest = ((st_i / m) * (M*M) + (end_j / m) * M + (st_k / m)) % num_threads; // % num_threads???
        MPI_Send(send_buf, message_len, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
    }
    if (end_k != N - 1) {
        int k = end_k - 1;
        int message_len = (end_j - st_j) * (end_i - st_i);
        double *send_buf = calloc(message_len, sizeof(double));
        for (int i = st_i; i < end_i; i++) {
            for (int j = st_j; j < end_j; j++) {
                send_buf[i*N + j] = A[i][j][k];
            }
        }
        int dest = ((st_i / m) * (M*M) + (st_j / m) * M + (end_k / m)) % num_threads; // % num_threads???
        MPI_Send(send_buf, message_len, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
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
    /* printf("Proc %d/%d: local_s = %f\n", rank, num_threads, local_s); */
}


int main(int argc, char **argv) {

    // ПАМЯТЬ УТЕКАЕТ!!!
    // ПАМЯТЬ УТЕКАЕТ!!!
    // ПАМЯТЬ УТЕКАЕТ!!!
    // ПАМЯТЬ УТЕКАЕТ!!!
    // ПАМЯТЬ УТЕКАЕТ!!!
    // ПАМЯТЬ УТЕКАЕТ!!!
    int status = MPI_Init(&argc, &argv);
    if (status) { 
        printf("MPI not supported\nError with code %d\n", status);
        MPI_Abort(MPI_COMM_WORLD, status);
        return status;
    };

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_threads);
    printf("i'm %d / %d\n", rank, num_threads);

    M = (N + m - 1) / m;

    double start, end;
    int st_i, st_j, st_k;
    
    if (rank == 0) {
        start = MPI_Wtime();
    }
    
    init();
    int it;
    for (it = 1; it <= itmax; it++)
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
        /* relax(); */
        /* printf("it=%4i   eps=%f\n", it, eps); */
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


