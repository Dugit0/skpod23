import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
from collections import defaultdict
import os
import pprint

prog_types = ['openmp_hyperplain', 'openmp_redblack']
mod_types = ['O0', 'O2', 'O3']
mods_dict = {j: {i: [] for i in prog_types} for j in mod_types}
for name in os.listdir('out'):
    if name.endswith('err'):
        continue
    with open(f'out/{name}') as f_inp:
        for line in f_inp:
            if line.strip().startswith('Time = '):
                # print(name, end="  :  ")
                # print(line.strip())
                time = float(line.split('=')[1])
                break
        else:
            if name.endswith("out"):
                time = float(15 * 60)
    for prog_type in prog_types:
        if name.startswith(prog_type) and name.endswith("out"):
            n, mod, proc = name.replace(prog_type, '').replace('.out', '').split('_')
            n = int(n)
            proc = int(proc)
            mods_dict[mod][prog_type].append((n, proc, time))

for mod in mod_types:
    for alg, orig_data in mods_dict[mod].items():
        print(mod, alg)
        sorted_data = sorted(orig_data, key=lambda a: (a[0], a[1])).copy()
        print(*sorted_data, sep='\n')
        values = sorted(list({i[0] for i in sorted_data}))
        num_proc = sorted(list({i[1] for i in sorted_data}))
        n = len(values)
        m = len(num_proc)
        print(f"{n = }, {m = }")
        matr = [[sorted_data[i*m + j][2] for j in range(m)] for i in range(n)]
        # data_frame = pd.DataFrame(matr, columns=num_proc, index=values)
        data_frame = pd.DataFrame(matr, columns=num_proc, index=values)
        # data_frame = pd.pivot_table(data_frame, values='V', index='Proc', columns='Vals')
        # data_frame.pivot('Proc', 'Vals')
        # print(data_frame)
        plt.title(alg)
        plt.figure(figsize=(15, 5))
        sns.heatmap(data_frame, cmap='Blues', annot=True, fmt='.2f')
        plt.xlabel('number of threads')
        plt.ylabel('N')
        plt.savefig(f'{alg}_{mod}.png')
        plt.close()



