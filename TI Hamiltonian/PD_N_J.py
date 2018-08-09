# This program is aimed to search the Topological Phase
# in 2D phase diagram of N(hybrdization gap) & J(Spin Gap)

import numpy as np
import os
import random
from scipy import linalg
import matplotlib as mpl
mpl.use("Agg") # 在服务器端没有DISPLAY资源，因此必须换成无交互的后端。默认为"TkAgg"。
from matplotlib import pyplot as plt


from json import dumps, loads
from itertools import product
from datetime import datetime
from functools import partial
from multiprocessing import Lock, Pool, Queue, Manager

from TI_Film import Hamiltonian, Eig
import TopoInvCalc as TIC

settings = {'num_lines': 31,
            'pos_tol':  1e-2,
            'gap_tol': 0.1,
            'move_tol': 0.1,
            'iterator': range(30, 51, 4),
            'min_neighbour_dist': 1e-4,
            }

def nt():
    return datetime.now().strftime("%y-%m-%d-%H-%M-%S")


def CalcGap():
    N_min, N_max, J = 6, 30, 0.00
    Gap = np.zeros([N_max+1], dtype=float)
    for N in range(N_min, N_max+1):
        h = Hamiltonian(N=N, J=J)   # By Default it is Spin-Z
        eig = np.array([x.real for x in (linalg.eig(h(0, 0))[0])])
        eig.sort()
        Gap[N] = eig[2*N]-eig[2*N-1]

    plt.subplot()
    plt.plot(list(range(N_min, N_max+1)), Gap[N_min:N_max+1]*1000)
    plt.xlabel(r"$N$")
    plt.ylabel(r"$Gap(\rm meV)$")
    plt.title("Gap vs Thickness, No Spin")
    plt.savefig("Gap_vs_Thickness_No_Spin.png")
    return


def TopoOrder(res, _Chern_tol=0.1):
    '''
    Accept an Result from TopoInvCalc,
    return
        0 if trivial;
        1 if Z2 but not Chern;
        2 if Chern but not Z2;
        3 if Chern and Z2;
    '''
    C, Z = 1 if abs(res.Chern) > _Chern_tol else 0, int(res.Z2)
    return C*2+Z


def Run(_N, _J, i, j, Phase):
    # 使用Manager().list在多进程通信
    # 在Python3.5下表现有问题：不报错，但是无法写入
    h = Hamiltonian(N=_N, J=_J)
    res = TIC.Calc(h, CalcZ2=True,settings=settings)
    print(res.Chern)
    Phase[i][j] = TopoOrder(res)
    return


def PhaseDiag():
    N_min, N_max, J_min, J_max, NJ = 6, 20, 0.00, 0.02, 20
    N = np.array(list(range(N_min, N_max+1)), dtype=int)
    J = np.linspace(J_min, J_max, num=NJ, endpoint=True)

    p, m = Pool(), Manager()
    Phase = m.list([m.list([0]*NJ) for i in range(N_max-N_min+1)])
    # 二维数组实在是太坑爹了......
    # 坑爹*2： 内层的m.list是对象，所以如果直接用[...]*N的方式写会导致写入引用

    for i, j in product(list(range(N_max-N_min+1)), list(range(NJ))):
        N_, J_ = N[i], J[j]
        p.apply_async(Run, args=(N_, J_, i, j, Phase,))

    p.close()
    p.join()
    P = [list(x) for x in list(Phase)]
    stime = nt()
    with open("PD_Data_"+stime+".txt", "w") as f:
        d = dict(N_min=N_min, N_max=N_max, J_min=J_min,
                 J_max=J_max, NJ=NJ, Phase=P)
        f.write(dumps(d))

    def Draw(filename=""):
        fig, ax = plt.subplots()
        cmap = mpl.colors.ListedColormap(["r", "g", "b", "c"])
        norm = mpl.colors.BoundaryNorm(list(range(5)), cmap.N)
        x, y = np.meshgrid(
            np.append(N, 2*N[-1]-N[-2]), np.append(J, 2*J[-1]-J[-2]))
        c = ax.pcolormesh(x, y, np.array(P).T, cmap=cmap, vmin=0, vmax=3)
        ax2 = fig.add_axes([0.92, 0.1, 0.03, 0.8])
        cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm)
        ax.set_title("Phase Diagram of N & J")
        if filename != "":
            plt.savefig(filename+".png")
        else:
            plt.show()
    Draw("PhaseDiagram_"+stime)


if __name__ == "__main__":
    # CalcGap()
    PhaseDiag()
