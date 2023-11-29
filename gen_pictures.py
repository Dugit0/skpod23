import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
from collections import defaultdict
import os
import pprint

prog_types = ['openmp_hyperplain', 'openmp_redblack']
data_dict = {i: [] for i in prog_types}
for name in os.listdir('out'):
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
            n, proc = map(int, name.replace(prog_type, '').replace('.out', '').split('_'))
            data_dict[prog_type].append((n, proc, time))

for alg, orig_data in data_dict.items():
    print(alg)
    sorted_data = sorted(orig_data, key=lambda a: (-a[0], a[1])).copy()
    # print(*sorted_data, sep='\n')
    values = sorted(list({i[0] for i in sorted_data}), reverse=True)
    num_proc = sorted(list({i[1] for i in sorted_data}))
    n = len(num_proc)
    m = len(values)
    matr = [[sorted_data[i*m + j][2] for j in range(m)] for i in range(n)]
    data_frame = pd.DataFrame(matr, columns=num_proc, index=values)
    # data_frame = pd.pivot_table(data_frame, values='V', index='Proc', columns='Vals')
    # data_frame.pivot('Proc', 'Vals')
    print(data_frame)
    plt.title(alg)
    sns.heatmap(data_frame, cmap='Blues', annot=True, fmt='.2f')
    plt.xlabel('number of threads')
    plt.ylabel('N')
    plt.savefig(f'{alg}.png')
    plt.close()
    # pprint.pprint(matr)



