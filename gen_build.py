import os
gen = 'gen_prog'
processers = [1, 2, 8, 20, 100]
with open("build.sh", "w") as f_out:
    f_out.write("#!/bin/bash\n")
    f_out.write(f"cd {gen}\n")
    for name in os.listdir(gen):
        f_out.write(f"xlc -qsmp=omp -O2 -o ../build_prog/{name.split('.')[0]} {name}\n")
os.chmod("build.sh", 0o744)
