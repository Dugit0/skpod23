#include <mpi.h>
#include <mpi-ext.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <time.h>

#define DEBUG 0

struct ErrorImitation {
    int flag;
    int rank;
    int iteration;
};

enum {
    BACKUP_SIZE = 2
};

ErrorImitation breakdowns[BACKUP_SIZE] = {{0, 2, 3}, {0, 1, 7}};
int break_number = 0;
int rank, size, hot_size;
MPI_Comm comm_world;

static void verbose_errhandler(MPI_Comm* comm, int* err, ...) {

    breakdowns[break_number].flag = 1;
    break_number++;
    int old_rank = rank;

    MPIX_Comm_revoke(comm_world);
    MPIX_Comm_shrink(comm_world, &comm_world);

    MPI_Comm_rank(comm_world, &rank);
    MPI_Comm_size(comm_world, &size);
    printf("%d/%d: Breakdown %d. old_rank = %d, new_rank = %d\n",
           rank, size, break_number,old_rank, rank);
    usleep(100);

    throw 0;
}

void my_func(int iter) {
    if (rank == 0) {
        int buf = 42;
        if (DEBUG >= 5)
            printf("%d/%d: Send 42 to %d.\n", rank, size, rank + 1);
        MPI_Send(&buf, 1, MPI_INT, rank + 1, 0, comm_world);
    } else if (rank == hot_size - 1) {
        int buf;
        MPI_Status status;
        MPI_Recv(&buf, 1, MPI_INT, rank - 1, 0, comm_world, &status);
        printf("i = %d, answer is %d\n", iter, buf);
    } else {
        int buf;
        MPI_Status status;
        MPI_Recv(&buf, 1, MPI_INT, rank - 1, 0, comm_world, &status);
        if (DEBUG >= 5)
            printf("%d/%d: Recv from %d, send to %d.\n", rank, size, rank - 1, rank + 1);
        MPI_Send(&buf, 1, MPI_INT, rank + 1, 0, comm_world);
    }
}

int main(int argc, char *argv[]) {
    MPI_Errhandler errh;
    MPI_Init(NULL, NULL);
    comm_world = MPI_COMM_WORLD;
    MPI_Comm_rank(comm_world, &rank);
    MPI_Comm_size(comm_world, &size);
    hot_size = size - BACKUP_SIZE;
    MPI_Comm_create_errhandler(verbose_errhandler, &errh);
    MPI_Comm_set_errhandler(comm_world, errh);

    printf("%d/%d: Created\n", rank, size);
    MPI_Barrier(comm_world);
    usleep(100);

    for (int i = 0; i < 10; i++) {
        try {
            if (!breakdowns[break_number].flag
                    && rank == breakdowns[break_number].rank
                    && i == breakdowns[break_number].iteration) {
                raise(SIGKILL);
            }
            if (rank < hot_size) {
                my_func(i);
            }
            usleep(100);
            MPI_Barrier(comm_world);
        } catch (int error) {
            i--;
            MPI_Barrier(comm_world);
        }
    }

    MPI_Barrier(comm_world);
    usleep(100);
    printf("%d/%d: Stayin' alive!\n", rank, size);
    MPI_Finalize();
}
