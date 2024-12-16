#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <time.h>


int broke_flag = 0;

static void verbose_errhandler(MPI_Comm* comm, int* err, ...) {
    int rank, size, len;
    char errstr[MPI_MAX_ERROR_STRING];
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Error_string( *err, errstr, &len );
    printf("Rank %d / %d: Notified of error %s\n",
           rank, size, errstr);
}

void my_func(int rank, int hot_size, int iter) {
    if (rank == 0) {
        int buf = 42;
        MPI_Send(&buf, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    } else if (rank == hot_size - 1) {
        int buf;
        MPI_Status status;
        MPI_Recv(&buf, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);
        printf("i = %d, answer is %d\n", iter, buf);
    } else {
        int buf;
        MPI_Status status;
        MPI_Recv(&buf, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);
        MPI_Send(&buf, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }
}

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Errhandler errh;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int hot_size = size - 1;
    MPI_Comm_create_errhandler(verbose_errhandler, &errh);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, errh);

    printf("%d/%d: Created\n", rank, size);
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100);

    for (int i = 0; i < 10; i++) {
        if (rank < hot_size) {
            my_func(rank, hot_size, i);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // if( rank == (size-1) ) raise(SIGKILL);
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100);
    printf("%d/%d: Stayin' alive!\n", rank, size);
    MPI_Finalize();
}
