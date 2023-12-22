import os
gen = 'gen_prog'
opt_flags = ['O0', 'O2', 'O3']
with open("build.sh", "w") as f_out:
    f_out.write("#!/bin/bash\n")
    f_out.write(f"cd {gen}\n")
    for name in os.listdir(gen):
        for opt_flag in opt_flags:
            # f_out.write(f"xlc -qsmp=omp -{opt_flag} -o ../build_prog/{name.split('.')[0]}_{opt_flag} {name}\n")
            f_out.write(f"gcc -fopenmp -std=c99 -{opt_flag} -o ../build_prog/{name.split('.')[0]}_{opt_flag} {name}\n")
os.chmod("build.sh", 0o744)
