import random
import time
import numpy as np
from functools import partial
from multiprocessing import Lock, Pool, Manager
from itertools import product

import matplotlib as mpl
mpl.use("Agg")
from matplotlib import pyplot as plt

random.seed()


def Run(i, Y, lock):
    print("Run: index=", i)
    time.sleep(int(random.random()*10))
    # lock.acquire()
    # print("Get Lock!, index=", i)
    Y[i] = i*i
    print("Write: index=", i)
    # lock.release()


def Run_without_lock(N_, J_, i, j, Z):
    # print("Run: index=", i, ",", j)
    # time.sleep(int(random.random()*10))
    # lock.acquire()
    # print("Get Lock!, index=", i)
    temp = int(random.random()*4)
    Z[i][j] = temp
    print("i,j,N,J,temp,Value=", i, j, N_, J_, temp, Z[i][j])
    # Z[i][j] = i*10+j

    # print("Write: index=", i, ",", j, " value=", Z[i][j])


def Draw(X, Y, Z):
    # print(X)
    # print(Y)
    # print(Z)
    x, y = np.meshgrid(X, Y, indexing="ij")

    fig, ax = plt.subplots()

    cmap = mpl.colors.ListedColormap(["r", "g", "b", "c"])
    norm = mpl.colors.BoundaryNorm(list(range(5)), cmap.N)
    c = ax.pcolormesh(x, y, np.array(Z), cmap=cmap, vmin=0, vmax=3)

    ax2 = fig.add_axes([0.92, 0.1, 0.03, 0.8])

    cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm)
    ax.set_title("Title")
    # fig.colorbar(cmap, ax=ax, norm=norm)
    # fig.colorbar(c, ax=ax)
    plt.savefig("TestForMP.png")
    # plt.show()
    return


if __name__ == "__main__":
    # p, m = Pool(), Manager()
    # # l = m.Lock()
    # # global Y
    # N_min, N_max, J_min, J_max, NJ = 15, 20, 0.00, 0.02, 5
    # N = np.array(list(range(N_min, N_max+1)), dtype=int)
    # J = np.linspace(J_min, J_max, num=NJ, endpoint=True)

    # x, y = np.meshgrid(N, J, indexing="ij")
    # Z = m.list([m.list([0 for j in range(NJ)])for i in range(N_max-N_min+1)])
    # # pRun = partial(Run, lock=l)
    # for i, j in product(list(range(N_max-N_min+1)), list(range(NJ))):
    #     p.apply_async(Run_without_lock, args=(N[i], J[j], i, j, Z,))
    # p.close()
    # p.join()
    # z = [list(x) for x in list(Z)]
    # with open("TestForMP.txt","w") as f:
    #     f.write(str(z))

    # # print(z)
    # Draw(N, J, z)
    A = np.linspace(1,3,4,endpoint=True)
    B = np.linspace(4,7,5,endpoint=True)
    C = [[1, 2, 3, 0], [0, 0, 3, 2], [2, 2, 1, 1]]
    Draw(A, B, C)
