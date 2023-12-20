/* #include "mpi.h" */
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
double check_sum = 0.;
double eps, local_eps;
int rank, num_threads;

double A[N][N][N];

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
    Slice inds;
    Triplet send_to;
    Triplet recv_from;
} Info;


void print_info(Info i) {
    printf("rank = %d, begin = (%d %d %d), end = (%d %d %d)\nsend_to = (%d %d %d), recv_from = (%d %d %d)",
            i.proc_rank,
            i.inds.beg_i, i.inds.beg_j, i.inds.beg_k,
            i.inds.end_i, i.inds.end_j, i.inds.end_k,
            i.send_to.i, i.send_to.j, i.send_to.k,
            i.recv_from.i, i.recv_from.j, i.recv_from.k);
}

/* int get_proc_from_cube(int cube_ind) { */
/* } */

Info *init_metadata() {
    Info *res = calloc(M*M*M, sizeof(*res));
    for (int i = 0; i < M*M*M; i++) {
        res[i].proc_rank = i % num_threads;
        res[i].inds.beg_i = i / (M*M) * m;
        res[i].inds.beg_j = i % (M*M) / M * m;
        res[i].inds.beg_k = i % M * m;
        res[i].inds.end_i = Min(res[i].inds.beg_i + m, N - 1);
        res[i].inds.end_j = Min(res[i].inds.beg_j + m, N - 1);
        res[i].inds.end_k = Min(res[i].inds.beg_k + m, N - 1);
        
        int buf;
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

    }
    
    return res;
}



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

void relax()
{
    return;
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

    M = (N + m - 1) / m;
    printf("N = %d, m = %d, M = %d\n", N, m, M);
    

    /* int status = MPI_Init(&argc, &argv); */
    /* if (status) {  */
    /*     printf("MPI not supported\nError with code %d\n", status); */
    /*     MPI_Abort(MPI_COMM_WORLD, status); */
    /*     return status; */
    /* }; */
    /*  */
    /* MPI_Comm_rank(MPI_COMM_WORLD, &rank); */
    /* MPI_Comm_size(MPI_COMM_WORLD, &num_threads); */
    /* if (debug >= 5) */
    /* printf("Created %d/%d\n", rank, num_threads); */
    /*  */
    /*  */
    /* double start, end; */
    /* int st_i, st_j, st_k; */
    /* MPI_Barrier(MPI_COMM_WORLD); */
    /*  */
    /* if (rank == 0) { */
    /*     start = MPI_Wtime(); */
    /* } */
    /*  */
    /* init(); */
    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* int it; */
    /* // for (it = 1; it <= itmax; it++) */
    /* for (it = 1; it <= 1; it++) */
    /* { */
    /*     eps = 0.; */
    /*     local_eps = 0.; */
    /*     for (int cur_cube = rank; cur_cube < M*M*M; cur_cube += num_threads) { */
    /*         st_i = cur_cube / (M*M) * m; */
    /*         st_j = (cur_cube % (M*M)) / M * m; */
    /*         st_k = (cur_cube % M) * m; */
    /*         relax(st_i, st_j, st_k, cur_cube / num_threads); */
    /*     } */
    /*     MPI_Allreduce(&local_eps, &eps, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); */
    /*     MPI_Barrier(MPI_COMM_WORLD); */
    /*     if (rank == 0) { */
    /*         printa(); */
    /*         printf("===============================\n"); */
    /*         printf("-------------------------------\n"); */
    /*         printf("===============================\n"); */
    /*     } */
    /*     if (eps < maxeps) { */
    /*         break; */
    /*     } */
    /* } */
    /*  */
    /* double s = 0.; */
    /* for (int cur_cube = rank; cur_cube < M*M*M; cur_cube += num_threads) { */
    /*     st_i = cur_cube / (M*M) * m; */
    /*     st_j = (cur_cube % (M*M)) / M * m; */
    /*     st_k = (cur_cube % M) * m; */
    /*     s += verify(st_i, st_j, st_k); */
    /* } */
    /* printf("Proc %d/%d: s = %f\n", rank, num_threads, s); */
    /* MPI_Reduce(&s, &check_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); */
    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* if (rank == 0) { */
    /*     end = MPI_Wtime(); */
    /*     printf("check_sum = %lf\n", check_sum); */
    /*     printf("Time = %f, num_threads = %d\n", end - start, num_threads); */
    /* } */
    /* MPI_Finalize(); */

    return 0;
}


