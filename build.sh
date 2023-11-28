cd gen_prog
xlc -qsmp=omp -O2 -o ../build_prog/openmp_for100 openmp_for100.c
xlc -qsmp=omp -O2 -o ../build_prog/openmp_for200 openmp_for200.c
xlc -qsmp=omp -O2 -o ../build_prog/openmp_for500 openmp_for500.c
xlc -qsmp=omp -O2 -o ../build_prog/openmp_for1000 openmp_for1000.c
xlc -qsmp=omp -O2 -o ../build_prog/openmp_redblack100 openmp_redblack100.c
xlc -qsmp=omp -O2 -o ../build_prog/openmp_redblack200 openmp_redblack200.c
xlc -qsmp=omp -O2 -o ../build_prog/openmp_redblack500 openmp_redblack500.c
xlc -qsmp=omp -O2 -o ../build_prog/openmp_redblack1000 openmp_redblack1000.c
xlc -qsmp=omp -O2 -o ../build_prog/prog_without_par100 prog_without_par100.c
xlc -qsmp=omp -O2 -o ../build_prog/prog_without_par200 prog_without_par200.c
xlc -qsmp=omp -O2 -o ../build_prog/prog_without_par500 prog_without_par500.c
xlc -qsmp=omp -O2 -o ../build_prog/prog_without_par1000 prog_without_par1000.c
xlc -qsmp=omp -O2 -o ../build_prog/var20100 var20100.c
xlc -qsmp=omp -O2 -o ../build_prog/var20200 var20200.c
xlc -qsmp=omp -O2 -o ../build_prog/var20500 var20500.c
xlc -qsmp=omp -O2 -o ../build_prog/var201000 var201000.c
