# This program is aimed to search the Topological Phase
# in 2D phase diagram of M_T & M_B

import numpy as np
import os
import random
from scipy import linalg
import matplotlib as mpl
mpl.use("Agg")  # 在服务器端没有DISPLAY资源，因此必须换成无交互的后端。默认为"TkAgg"。
from matplotlib import pyplot as plt

# from TI_AFM import Ham_Small as Hamiltonian
from TI_AFM import Hamiltonian

from json import dumps, loads
from itertools import product
from datetime import datetime
from functools import partial
from multiprocessing import Lock, Pool, Queue, Manager
# 使用Manager().list在多进程通信
# 在Python3.5下表现有问题：不报错，但是无法写入
import TopoInvCalc as TIC
from PD_N_J import nt


settings = {'num_lines': 51,
            'pos_tol':  1e-2,
            'gap_tol': 0.05,
            'move_tol': 0.2,
            'iterator': range(50, 81, 4),
            'min_neighbour_dist': 1e-4,
            }


def Draw(T, B, P, title="Phase Diagram of M_T & M_B", filename=""):
    fig, ax = plt.subplots()
    cmap = mpl.colors.ListedColormap(["r", "g", "b", "c"])
    norm = mpl.colors.BoundaryNorm(list(range(5)), cmap.N)
    # x, y = np.meshgrid(
    # np.append(T, 2*T[-1]-T[-2]), np.append(B, 2*B[-1]-B[-2]))
    x, y = np.meshgrid(T, B)
    c = ax.pcolormesh(x, y, np.array(P).T, cmap=cmap, vmin=0, vmax=3)
    ax2 = fig.add_axes([0.92, 0.1, 0.03, 0.8])
    cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm)
    ax.set_title(title)
    if filename != "":
        plt.savefig(filename+".png")
    else:
        plt.show()


def Run_1(M_T, M_B, i, j, Phase):
    '''
    This func calculate the system with spin distribution: +z in the upper half and -z in the lower half.
    If the # of layers is odd, then the layer in the very middle will have no spin to make sure the total spin in Z direction of this system to be ZERO.
    '''
    print("Start Calculation: M_T=%0.3f , M_B=%.3f" % (M_T, M_B))

    h = Hamiltonian(M_T, M_B)
    res = TIC.Calc(h, CalcZ2=True, LogOut=False, settings=settings)
    Phase[i][j] = TIC.TopoOrder(res, _Chern_tol=0.3)
    print("End Calculation: M_T=%0.3f , M_B=%.3f, Result: C=%.4f , Z2=%s" %
          (M_T, M_B, res.Chern, str(res._Z2)))
    return


def PhaseDiag(func, title="Phase Diagram of M_T & M_B(Large Hamiltonian)"):
    T_start = datetime.now()
    print("Start Calculation at ", str(T_start))

    TRange, BRange, N_T, N_B = 0.2, 0.2, 10, 10
    T = np.linspace(-TRange, TRange, num=N_T, endpoint=True)
    B = np.linspace(-BRange, BRange, num=N_B, endpoint=True)

    p, m = Pool(), Manager()
    Phase = m.list([m.list([0]*(N_B)) for i in range(N_T)])
    # 二维数组实在是太坑爹了......
    # 坑爹*2： 内层的m.list是对象，所以如果直接用[...]*N的方式写会导致写入引用

    for i, j in product(list(range(N_T)), list(range(N_B))):
        M_T, M_B = T[i], B[j]
        p.apply_async(func, args=(M_T, M_B, i, j, Phase,))

    p.close()
    p.join()
    P = [list(x) for x in list(Phase)]
    stime = nt()
    with open("PD_TB_Data_"+stime+".txt", "w") as f:
        d = dict(TRange=TRange, BRange=BRange,
                 N_T=N_T, N_B=N_B, Phase=P)
        f.write(dumps(d))

    Draw(T, B, P, title=title,
         filename="PhaseDiag_TB_"+stime)
    T_end = datetime.now()
    print("End Calculation at ", str(T_end))
    print("Total time:", str(T_end-T_start))


if __name__ == "__main__":
    # CalcGap()
    PhaseDiag(Run_1)
