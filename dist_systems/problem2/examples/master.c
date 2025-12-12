#include <stdio.h>
#include "mpi.h"
int main(int argc, char **argv)
{
 int size, rank1, rank2;
 MPI_Status status;
 MPI_Comm intercomm;
 char slave[10]="./slave";
 MPI_Init(&argc, &argv);
 MPI_Comm_size(MPI_COMM_WORLD, &size);
 MPI_Comm_spawn(slave, MPI_ARGV_NULL, 2, MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomm, MPI_ERRCODES_IGNORE);
 MPI_Recv(&rank1, 1, MPI_INT, 0, 0, intercomm, &status);
 MPI_Recv(&rank2, 1, MPI_INT, 1, 1, intercomm, &status);
 printf("Slaves %d and %d are working\n", rank1, rank2);
 MPI_Finalize();
 return 0;
}
