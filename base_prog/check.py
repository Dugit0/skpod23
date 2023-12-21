with open('mpi_out.txt') as f_inp:
    mpi = f_inp.read()
with open('real_out.txt') as f_inp:
    real = f_inp.read()

def prepare(s):
    s = sorted(list(map(float, s.strip().split('\n'))))
    return s

mpi = prepare(mpi)
real = prepare(real)

print(len(mpi), len(real))
for i, j in zip(mpi, real):
    print(i, j)
