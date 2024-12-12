#!/usr/bin/env bash


function cmcmpi() {
    outfile=${1%%.*}
    echo ${outfile}
    mpic++ -Wall -Wextra -fsanitize=address,undefined,signed-integer-overflow,pointer-compare,pointer-subtract,leak,bounds,pointer-overflow -o "${outfile}.out" ${1}
}

echo "======== COMPILE ========" &&
mpic++ -Wall -Wextra -std=c++23 -o main.out main.cpp &&
# cmcmpi main.c &&
echo "========== RUN ==========" &&
mpirun -np 5 --use-hwthread-cpus ./main.out
# mpirun -np 5 ./main.out

# echo "============ ============"
