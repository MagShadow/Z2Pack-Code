# This program is aimed to search the Topo Phase in 2D
# phase diagram of J(Spin Exchange Gap) & Angle(Spin
# direction), basicly have 3 spin alignment:
#   1. Mirror Sym
#   2. Inversion Sym
#   3. Symmetric under combination of Mirror & Inversion

import numpy as np
import os
import random
from scipy import linalg
import matplotlib as mpl
mpl.use("Agg")  # 在服务器端没有DISPLAY资源，因此必须换成无交互的后端。默认为"TkAgg"。
from matplotlib import pyplot as plt

from json import dumps, loads
from itertools import product
from datetime import datetime
from functools import partial
from multiprocessing import Lock, Pool, Queue, Manager
# 使用Manager().list在多进程通信
# 在Python3.5下表现有问题：不报错，但是无法写入
from TI_Film import Hamiltonian, Eig
from PD_N_J import nt
import TopoInvCalc as TIC

settings = {'num_lines': 51,
            'pos_tol':  1e-2,
            'gap_tol': 0.05,
            'move_tol': 0.2,
            'iterator': range(50, 81, 4),
            'min_neighbour_dist': 1e-4,
            }
N = 12


def Draw(T, J, P, title="Phase Diagram of N & J", filename=""):
    fig, ax = plt.subplots()
    cmap = mpl.colors.ListedColormap(["r", "g", "b", "c"])
    norm = mpl.colors.BoundaryNorm(list(range(5)), cmap.N)
    # x, y = np.meshgrid(
    #     np.append(N, 2*N[-1]-N[-2]), np.append(J, 2*J[-1]-J[-2]))
    x, y = np.meshgrid(T, J)
    c = ax.pcolormesh(x, y, np.array(P).T, cmap=cmap, vmin=0, vmax=3)
    ax2 = fig.add_axes([0.92, 0.1, 0.03, 0.8])
    cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm)
    ax.set_title(title)
    if filename != "":
        plt.savefig(filename+".png")
    else:
        plt.show()


def Run_1(_T, _J, i, j, Phase, n=None):
    '''
    This func calculate the system with Mirror Symmetry,   which means that for the upper half and the lower half,the parallel components changes sign while 
    :para `n` defines how many layers (in upper half) contains the Spin term. In default, all of the upper half have the spin term.
    '''
    print("Start Calculation: Theta=%.3f , J=%.3f" % (_T, _J))
    _S = np.zeros([N, 3])
    layer = int(N/2) if n == None else min(int(n), int(N/2))
    for i in range(layer):
        _S[i, 0], _S[N-i-1, 0] = 1, 1
        _S[i, 1], _S[N-i-1, 1] = _T, _T
        _S[N-i-1, 2] = np.pi
    S = np.array([([s[0]*np.sin(s[1])*np.cos(s[2]), s[0] *
                    np.sin(s[1]) * np.sin(s[2]), s[0]*np.cos(s[1])])for s in _S])
    print(S)

    h = Hamiltonian(N=N, J=_J, S=S)
    res = TIC.Calc(h, CalcZ2=True, LogOut=False, settings=settings)
    Phase[i][j] = TIC.TopoOrder(res)
    print("End Calculation: N=%d , J=%.3f, Result: C=%.4f , Z2=%s" %
          (N, _J, res.Chern, str(res._Z2)))
    return


def PhaseDiag(func, title="Phase Diagram of Theta & J"):
    T_start = datetime.now()
    print("Start Calculation at ", str(T_start))

    J_min, J_max, N_J = 0.00, 0.04, 11
    Theta_min, Theta_max, N_Theta = 0, np.pi, 11
    J = np.linspace(J_min, J_max, N_J, endpoint=True)
    T = np.linspace(Theta_min, Theta_max, N_Theta, endpoint=True)

    # Run Calculation in Multi Process
    p, m = Pool(), Manager()
    Phase = m.list([m.list([0]*N_J) for i in range(N_Theta)])

    for i, j in product(list(range(N_Theta)), list(range(N_J))):
        T_, J_ = T[i], J[j]
        p.apply_async(func, args=(T_, J_, i, j, Phase,))

    p.close()
    p.join()

    P = [list(x) for x in list(Phase)]
    stime = nt()
    with open("PD_Data_"+stime+".txt", "w") as f:
        d = dict(Theta_min=Theta_min, Theta_max=Theta_max, N_Theta=N_Theta, J_min=J_min,
                 J_max=J_max, N_J=N_J, Phase=P)
        f.write(dumps(d))

    Draw(T, J, P, title=title,
         filename="PhaseDiag_"+stime)
    T_end = datetime.now()
    print("End Calculation at ", str(T_end))
    print("Total time:", str(T_end-T_start))


if __name__ == "__main__":
    PhaseDiag(Run_1,title="PhaseDiag: Mirror Symmetry")
