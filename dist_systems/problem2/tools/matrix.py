import numpy as np

def vizualize(m):
    print(m)

N = 6 
"""
https://stackoverflow.com/questions/53173098/diagonal-snake-filling-array
"""
a = np.zeros((N, N))

from itertools import count, islice
def my_gen(N): # generator of terms
    for n in range(0, 2*N - 1):
        if n >= N:
            n = n - N
            for m in range(n+1):
                yield (N - m - 1, N- (n-m) - 1) if n % 2 else (N - (n-m) - 1, N - m - 1)
        else:
            for m in range(n+1):
                yield (m, n-m) if n % 2 else (n-m, m)

def triple_gen(N):
    for k in range(0, N):
        for ki in range(k + 1):
            for a in range(ki + 1):
                print((a, ki-a, k - ki))
        print()
    print("-----------")
    for k in range(N, 2*N):
        k = k - N
        for ki in range(k + 1):
            for a in range(ki + 1):
                print((N - a - 1, N - (ki-a) - 1, k - ki))
        print()

triple_gen(N)

cnt = 1
gen = my_gen(N)
for i in range(N*N):
    try:
        tmp = next(gen)
        #print(tmp)
    except StopIteration:
        pass
    a[tmp] = cnt
    cnt += 1
print(a)

