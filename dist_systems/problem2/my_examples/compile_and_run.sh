#!/usr/bin/env bash
set -e

# function cmcmpi() {
#     outfile=${1%%.*}
#     echo ${outfile}
#     mpic++ -Wall -Wextra -fsanitize=address,undefined,signed-integer-overflow,pointer-compare,pointer-subtract,leak,bounds,pointer-overflow -o "${outfile}.out" ${1}
# }

echo "====== COMPILE ======"
mpic++ -Wall -Wextra -std=c++23 -o example.out example.cpp
echo "======== RUN ========"
mpirun -np 8 --use-hwthread-cpus --with-ft ulfm ./example.out

# echo "============ ============"
