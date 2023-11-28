import os
base = 'base_prog'
gen = 'gen_prog'
build = 'build_prog'
conf = 'configs'
processers = [1, 2, 8, 20, 100]
home = '/home_edu/edu-cmc-skpod23-321/edu-cmc-skpod23-321-01/'

for name in os.listdir(build):
    for proc in processers:
        with open(f"configs/{name}_{proc}", "w") as f:
            f.write("#BSUB -W 00:15\n")
            f.write(f"#BSUB -o {home}out/{name}_{proc}.out\n")
            f.write(f"#BSUB -e {home}out/{name}_{proc}.err\n")
            f.write('\n')
            f.write("#BSUB -n {}\n".format(proc//8 + 1))
            f.write(f"OMP_NUM_THREADS={proc} {name}\n")


