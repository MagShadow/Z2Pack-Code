# This program is aimed to search the Topological Phase
# in 2D phase diagram of N(hybrdization gap) & J(Spin Gap)

import numpy as np
import os
import random
from scipy import linalg
from matplotlib import pyplot as plt
import matplotlib as mpl

from json import dumps, loads
from itertools import product
from datetime import datetime
from functools import partial
from multiprocessing import Lock, Pool, Queue, Manager

from TI_Film import Hamiltonian, Eig
import TopoInvCalc as TIC


def nt():
    return datetime.now().strftime("%y-%m-%d-%H-%M-%S")


def CalcGap():
    N_min, N_max, J = 6, 30, 0.00
    Gap = np.zeros([N_max+1], dtype=float)
    for N in range(N_min, N_max+1):
        # S_ = np.array([[1, 0, 0]]*N)
        # S = np.array([([s[0]*np.sin(s[1])*np.cos(s[2]), s[0]*np.sin(s[1])
        #                 * np.sin(s[2]), s[0]*np.cos(s[1])])for s in S_])
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


# def init_lock(l):
#     global lock
#     lock = l
def Run(_N, _J, i, j, Phase):
    h = Hamiltonian(N=_N, J=_J)
    res = TIC.Calc(h, CalcZ2=True)
    Phase[i][j] = TopoOrder(res)
    # Phase[i][j] = int(random.random()*4)
    # print("New Run!")
    # lock.aquire()
    # print("Get Lock"+datetime.now().strftime("%y-%m-%d-%H-%M-%S"))
    # Phase[i][j] = i+j
    # print("Write:(", i, ",", j, ")=", Phase[i][j])


def PhaseDiag():
    N_min, N_max, J_min, J_max, NJ = 15, 25, 0.00, 0.02, 10
    N = np.array(list(range(N_min, N_max+1)), dtype=int)
    # for n in N:
    #     print(type(n))
    J = np.linspace(J_min, J_max, num=NJ, endpoint=True)
    x, y = np.meshgrid(N, J, indexing="ij")
    # global Phase
    # Phase = np.zeros([N_max+1, NJ], dtype=int)
    # lock = Lock()

    p, m = Pool(), Manager()
    Phase = m.list([m.list([0]*NJ) for i in range(N_max-N_min+1)])
    # 二维数组实在是太坑爹了......
    # 坑爹*2

    # P = [list(x) for x in list(Phase)]
    # P = [list(x) for x in list(Phase)[N_min:N_max+1]]
    # print(P)
    # print()
    # print(Phase)
    # l = m.Lock()

    # print("Release Lock"+datetime.now().strftime("%y-%m-%d-%H-%M-%S"))
    # lock.release()

    # p_Run = partial(Run, lock=l)
    for i, j in product(list(range(N_max-N_min+1)), list(range(NJ))):
        # print(i,j,N[i],J[j])
        # In Default Settings, Spin is in z-axis
        # h = Hamiltonian(N=N[i], J=J[j])
        N_, J_ = N[i], J[j]
        p.apply_async(Run, args=(N_, J_, i, j, Phase,))

    p.close()
    p.join()
    # Res=[[]]
    P = [list(x) for x in list(Phase)]
    # print(P)
    # pp = [list(x) for x in list(Phase)]
    # print(pp)
    stime = nt()
    with open("PD_Data_"+stime+".txt", "w") as f:
        d = dict(N_min=N_min, N_max=N_max, J_min=J_min,
                 J_max=J_max, NJ=NJ, Phase=P)
        f.write(dumps(d))

    fig, ax = plt.subplots()
    cmap = mpl.colors.ListedColormap(["r", "g", "b", "c"])
    norm = mpl.colors.BoundaryNorm(list(range(5)), cmap.N)
    x, y = np.meshgrid(
        np.append(N, 2*N[-1]-N[-2]), np.append(J, 2*J[-1]-J[-2]))
    # print(np.append(J, 2*J[-1]-J[-2]))
    c = ax.pcolormesh(x, y, np.array(P).T, cmap=cmap, vmin=0, vmax=3)
    ax2 = fig.add_axes([0.92, 0.1, 0.03, 0.8])
    cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm)
    ax.set_title("Phase Diagram of N & J")
    plt.savefig("PhaseDiag_"+stime+".png")
    # plt.show()


if __name__ == "__main__":
    # CalcGap()
    PhaseDiag()

    # print(Phase)
