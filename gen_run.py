import os
with open("run.sh", "w") as f_out:
    f_out.write(f'#!/bin/bash\n')
    for name in os.listdir('configs'):
        f_out.write(f'bsub < configs/{name}\n')
os.chmod("run.sh", 0o744)
