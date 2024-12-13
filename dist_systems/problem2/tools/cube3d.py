import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from itertools import cycle, product, permutations
from tqdm import tqdm

n = 50
# Create axis
axes = [n, n, n]

# Create Data
data = np.ones(axes)

# control Transparency
alpha = 0.7

# control colour
colors = np.empty(axes + [4])

colors[:] = [1, 1, 1, alpha] # red


print("Create c:")
c = list()
for i in tqdm(product("01", repeat=3)):
    c.append([int(j) for j in i])
c.pop()
c = cycle(c)


def triple_gen(N):
    sqrt_from_2_approx = 1.4
    #len of cube diagonal == sqrt(2)*N
    #we take resized cube, and fill lower
    #triangle for our resized cube =>
    #our one-sized cube will be filled fully
    #for k in range(0, int(2*sqrt_from_2_approx*N)):
    resized_length = 3 * N - 2
    print("Triple gen")
    for k in tqdm(range(0, resized_length)):
        col = next(c)
        for ki in range(k + 1):
            for a in range(ki + 1):
                if a < N and (ki-a) < N and (k - ki) < N:
                    #do stuff
                    colors[a, ki-a, k - ki] = col + [alpha]


triple_gen(n)

# Plot figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Voxels is used to customizations of the
# sizes, positions and colors.
ax.voxels(data, facecolors=colors, edgecolors='grey')
plt.savefig("out.png", dpi=300)
