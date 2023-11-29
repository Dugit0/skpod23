# import os
# build = 'build_prog'
# threads = [1, 2, 8, 20, 100]
# with open("build.sh", "w") as f_out:
#     f_out.write("#!/bin/bash\n")
#     f_out.write(f"cd {gen}\n")
#     for name in os.listdir(gen):
#         f_out.write(f"xlc -qsmp=omp -O2 -o ../build_prog/{name.split('.')[0]} {name}\n")
import os
with open("run.sh", "w") as f_out:
    f_out.write(f'#!/bin/bash\n')
    for name in os.listdir('config'):
        f_out.write(f'bsub < config/{name}\n')
