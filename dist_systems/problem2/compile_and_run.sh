#!/usr/bin/env bash
set -e

# function cmcmpi() {
#     outfile=${1%%.*}
#     echo ${outfile}
#     mpic++ -Wall -Wextra -fsanitize=address,undefined,signed-integer-overflow,pointer-compare,pointer-subtract,leak,bounds,pointer-overflow -o "${outfile}.out" ${1}
# }

echo "====== COMPILE dist_syst_mpi.cpp ======"
mpic++ -Wall -Wextra -std=c++23 -o dist_syst_mpi.out dist_syst_mpi.cpp
echo "======== RUN dist_syst_mpi.cpp ========"
mpirun -np 10 --use-hwthread-cpus --with-ft ulfm ./dist_syst_mpi.out

# echo "============ ============"
