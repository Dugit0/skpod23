#include "mpi.h"
/* #include <math.h> */
#include <stdlib.h>
#include <stdio.h>
#define Max(a, b) (((a) > (b)) ? (a) : (b))
#define Min(a, b) (((a) < (b)) ? (a) : (b))
#define  debug 10


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
    Triplet send_to;
    Triplet recv_from;
} Info;


int N = 10;
int M;
int m = 5;
double maxeps = 0.1e-7;
int itmax = 100;
double check_sum = 0.;
double eps, local_eps;
int rank, num_threads;
Info *metadata;


Info *init_metadata() {
    Info *res = calloc(M*M*M, sizeof(*res));
    for (int i = 0; i < M*M*M; i++) {
        res[i].proc_rank = i % num_threads;
        res[i].beg.i = i / (M*M) * m;
        res[i].beg.j = i % (M*M) / M * m;
        res[i].beg.k = i % M * m;
        res[i].end.i = Min(res[i].beg.i + m, N - 1);
        res[i].end.j = Min(res[i].beg.j + m, N - 1);
        res[i].end.k = Min(res[i].beg.k + m, N - 1);
        
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



double ***init(Triplet beg, Triplet end) {
    int len_i = end.i - beg.i + 1;
    int len_j = end.j - beg.j + 1;
    int len_k = end.k - beg.k + 1;
    
    double ***res = calloc(len_i, sizeof(*res));

    for (int i = 0; i < len_i; i++) {
        res[i] = calloc(len_j, sizeof(*(res[i])));
        for (int j = 0; j < len_j; j++) {
            res[i][j] = calloc(len_k, sizeof(*(res[i][j])));
            for (int k = 0; k < len_k; k++) {
                int real_i = beg.i + i;
                int real_j = beg.j + j;
                int real_k = beg.k + k;
                if (real_i == 0 || real_i == N - 1 || real_j == 0 || real_j == N - 1 || real_k == 0 || real_k == N - 1) {
                    res[i][j][k] = 0.;
                } else {
                    res[i][j][k] = (4. + real_i + real_j + real_k);
                }
            }
        }
    }
    return res;
}

void free_cube(Triplet beg, Triplet end, double ***cube) {

}

void relax(Triplet beg, Triplet end, double ***cube)
{
    return;
}


double verify(Triplet beg, Triplet end, double ***cube) {
    int len_i = end.i - beg.i + 1;
    int len_j = end.j - beg.j + 1;
    int len_k = end.k - beg.k + 1;
    double local_s = 0.;
    for (int i = 0; i < len_i; i++) {
        for (int j = 0; j < len_j; j++) {
            for (int k = 0; k < len_k; k++) {
                int real_i = beg.i + i;
                int real_j = beg.j + j;
                int real_k = beg.k + k;
                local_s += cube[i][j][k] * (real_i + 1) * (real_j + 1) * (real_k + 1) / (N * N * N);
            }
        }
    }
    return local_s;
}


void print_cube(Triplet beg, Triplet end, double ***cube) {
    int len_i = end.i - beg.i + 1;
    int len_j = end.j - beg.j + 1;
    int len_k = end.k - beg.k + 1;
    for (int i = 0; i < len_i; i++) {
        for (int j = 0; j < len_j; j++) {
            for (int k = 0; k < len_k; k++) {
                printf("%f ", cube[i][j][k]);
            }
            printf("\n");
        }
        printf("=======================\n");
    }
}


void print_info(Info i) {
    printf("rank = %d, begin = (%d %d %d), end = (%d %d %d)\nsend_to = (%d %d %d), recv_from = (%d %d %d)",
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
            /* printf("\n--------------\n"); */
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



int main(int argc, char **argv) {

    m = 4;
    M = (N + m - 1) / m;
    num_threads = 8;
    printf("N = %d, m = %d, M = %d\n", N, m, M);
    Info *metadata = init_metadata();
    print_metadata(metadata);

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
    
    
    double start, end;
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        start = MPI_Wtime();
    }

    int len_cubes = (M*M*M - rank) / num_threads + 1;
    int *nums_cubes = calloc(len_cubes, sizeof(*nums_cubes));
    double ****cubes = calloc(len_cubes, sizeof(*cubes));

    int tmp_ind;
    tmp_ind = 0;
    for (int cur_cube = rank; cur_cube < M*M*M; cur_cube += num_threads) {
        Info cur_metadata = metadata[cur_cube];
        nums_cubes[tmp_ind] = cur_cube;
        cubes[tmp_ind] = init(cur_metadata.beg, cur_metadata.end);
    }
    MPI_Barrier(MPI_COMM_WORLD);
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
    if (rank == 0) {
        end = MPI_Wtime();
        printf("check_sum = %lf\n", check_sum);
        printf("Time = %f, num_threads = %d\n", end - start, num_threads);
    }
    MPI_Finalize();
    free(metadata);

    return 0;
}


