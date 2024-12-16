#!/usr/bin/env bash
set -e

echo "=== COMPILE prog_without_par.c ==="
gcc -Wall -Wextra -O2 -std=c18 -fopenmp -lm -o prog_without_par.out prog_without_par.c
echo "===== RUN prog_without_par.c ====="
./prog_without_par.out

echo "====== COMPILE skpod_mpi.c ======"
mpicc -Wall -Wextra -std=c18 -o skpod_mpi.out skpod_mpi.c
echo "======== RUN skpod_mpi.c ========"
mpirun -np 8 --use-hwthread-cpus ./skpod_mpi.out

echo "====== COMPILE dist_syst_mpi.cpp ======"
mpic++ -Wall -Wextra -std=c++23 -o dist_syst_mpi.out dist_syst_mpi.cpp
echo "======== RUN dist_syst_mpi.cpp ========"
mpirun -np 10 --use-hwthread-cpus ./dist_syst_mpi.out


# echo "============ ============"
# mpirun -n 8 valgrind --suppressions=/usr/local/share/openmpi/openmpi-valgrind.supp ./skpod_mpi.out
