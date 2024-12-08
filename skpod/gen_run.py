import os
completed_tests = os.listdir('out')
completed_tests = set(map(lambda a: a[:-4], completed_tests))
with open("run.sh", "w") as f_out:
    f_out.write(f'#!/bin/bash\n')
    for name in os.listdir('configs'):
        if name not in completed_tests:
            f_out.write(f'bsub < configs/{name}\n')
os.chmod("run.sh", 0o744)
