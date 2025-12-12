#include <stdio.h>
#include "mpi.h"
int main(int argc, char **argv)
{
 int rank, size;
 MPI_Comm intercomm;
 MPI_Init(&argc, &argv);
 MPI_Comm_get_parent(&intercomm);
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 MPI_Send(&rank, 1, MPI_INT, 0, rank, intercomm);
 MPI_Finalize();
 return 0;
}