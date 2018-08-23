# This program is aimed to search the Topological Phase
# in 2D phase diagram of N(hybrdization gap) & J(Spin Gap)

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
import TopoInvCalc as TIC

settings = {'num_lines': 51,
            'pos_tol':  1e-2,
            'gap_tol': 0.05,
            'move_tol': 0.2,
            'iterator': range(50, 81, 4),
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


def Draw(N, J, P, title="Phase Diagram of N & J", filename=""):
    fig, ax = plt.subplots()
    cmap = mpl.colors.ListedColormap(["r", "g", "b", "c"])
    norm = mpl.colors.BoundaryNorm(list(range(5)), cmap.N)
    x, y = np.meshgrid(
        np.append(N, 2*N[-1]-N[-2]), np.append(J, 2*J[-1]-J[-2]))
    c = ax.pcolormesh(x, y, np.array(P).T, cmap=cmap, vmin=0, vmax=3)
    ax2 = fig.add_axes([0.92, 0.1, 0.03, 0.8])
    cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm)
    ax.set_title(title)
    if filename != "":
        plt.savefig(filename+".png")
    else:
        plt.show()


def Run_0(_N, _J, _i, _j, Phase):
    '''
    This func calculate the system with uniform spin distribution in Z-direction.
    '''
    h = Hamiltonian(N=_N, J=_J)
    res = TIC.Calc(h, CalcZ2=True, LogOut=False, settings=settings)
    print(res.Chern)
    Phase[_i][_j] = TIC.TopoOrder(res)
    return


def Run_1(_N, _J, _i, _j, Phase):
    '''
    This func calculate the system with spin distribution: +z in the upper half and -z in the lower half.
    If the # of layers is odd, then the layer in the very middle will have no spin to make sure the total spin in Z direction of this system to be ZERO.
    '''
    print("Start Calculation: N=%d , J=%.3f" % (_N, _J))
    _S = np.zeros([_N, 3])
    for i in range(int(_N/2)):
        _S[i, 0], _S[_N-i-1, 0] = 1, -1
    S = np.array([([s[0]*np.sin(s[1])*np.cos(s[2]), s[0] *
                    np.sin(s[1]) * np.sin(s[2]), s[0]*np.cos(s[1])])for s in _S])
    h = Hamiltonian(N=_N, J=_J, S=S)
    res = TIC.Calc(h, CalcZ2=True, LogOut=False, settings=settings)
    Phase[_i][_j] = TIC.TopoOrder(res)
    print("End Calculation: N=%d , J=%.3f, Result: C=%.4f , Z2=%s" %
          (_N, _J, res.Chern, str(res._Z2)))
    # Phase[_i][_j] = _i*10+_j

    return


def Run_2(_N, _J, _i, _j, Phase):
    '''
    This func calculate the system with spin distribution: +x in the upper half and -x in the lower half.
    If the # of layers is odd, then the layer in the very middle will have no spin to make sure the total spin in Z direction of this system to be ZERO.
    '''
    print("Start Calculation: N=%d , J=%.3f" % (_N, _J))
    _S = np.zeros([_N, 3])
    for i in range(int(_N/2)):
        _S[i, 0], _S[i, 1] = 1, np.pi/2
        _S[_N-i-1, 0], _S[_N-i-1, 1] = 1, np.pi/2
        _S[i, 2], _S[_N-i-1, 2] = 0, np.pi
    S = np.array([([s[0]*np.sin(s[1])*np.cos(s[2]), s[0] *
                    np.sin(s[1]) * np.sin(s[2]), s[0]*np.cos(s[1])])for s in _S])
    # print(S)
    h = Hamiltonian(N=_N, J=_J, S=S)
    res = TIC.Calc(h, CalcZ2=True, LogOut=False, settings=settings)
    Phase[_i][_j] = TIC.TopoOrder(res)
    print("End Calculation: N=%d , J=%.3f, Result: C=%.4f , Z2=%s" %
          (_N, _J, res.Chern, str(res._Z2)))
    return


def Run_3(_N, _J, _i, _j, Phase, n=1):
    r'''
    This func calculate the system with spin distribution: +z in the first $n$ layer and -z in the last $n$ layer.

    '''
    print("Start Calculation: N=%d , J=%.3f" % (_N, _J))
    _S = np.zeros([_N, 3])
    assert n <= int(_N/2), "n should be less than half of the layer #。"
    for i in range(n):
        _S[i, 0], _S[_N-i-1, 0] = 1, -1
    S = np.array([([s[0]*np.sin(s[1])*np.cos(s[2]), s[0] *
                    np.sin(s[1]) * np.sin(s[2]), s[0]*np.cos(s[1])])for s in _S])
    # print(S)
    h = Hamiltonian(N=_N, J=_J, S=S)
    res = TIC.Calc(h, CalcZ2=True, LogOut=False, settings=settings)
    Phase[_i][_j] = TIC.TopoOrder(res)
    print("End Calculation: N=%d , J=%.3f, Result: C=%.4f , Z2=%s" %
          (_N, _J, res.Chern, str(res._Z2)))
    return


def Run_4(_N, _J, _i, _j, Phase, n=1):
    r'''
    This func calculate the system with spin distribution: +z in the first $n$ layer and the last $n$ layer.

    '''
    print("Start Calculation: N=%d , J=%.3f" % (_N, _J))
    _S = np.zeros([_N, 3])
    assert n <= int(_N/2), "n should be less than half of the layer #。"
    for i in range(n):
        _S[i, 0], _S[_N-i-1, 0] = 1, 1
    S = np.array([([s[0]*np.sin(s[1])*np.cos(s[2]), s[0] *
                    np.sin(s[1]) * np.sin(s[2]), s[0]*np.cos(s[1])])for s in _S])
    # print(S)
    h = Hamiltonian(N=_N, J=_J, S=S)
    res = TIC.Calc(h, CalcZ2=True, LogOut=False, settings=settings)
    Phase[_i][_j] = TIC.TopoOrder(res)
    print("End Calculation: N=%d , J=%.3f, Result: C=%.4f , Z2=%s" %
          (_N, _J, res.Chern, str(res._Z2)))
    return


def PhaseDiag(func, title="Phase Diagram of N & J"):
    T_start = datetime.now()
    print("Start Calculation at ", str(T_start))
    N_min, N_max, J_min, J_max, NJ = 6, 20, 0.00, 0.04, 20
    # N_min, N_max, J_min, J_max, NJ = 19, 21, 0.00, 0.02, 3
    N = np.array(list(range(N_min, N_max+1)), dtype=int)
    J = np.linspace(J_min, J_max, num=NJ, endpoint=True)

    p, m = Pool(), Manager()
    Phase = m.list([m.list([0]*NJ) for i in range(N_max-N_min+1)])
    # 二维数组实在是太坑爹了......
    # 坑爹*2： 内层的m.list是对象，所以如果直接用[...]*N的方式写会导致写入引用

    for i, j in product(list(range(N_max-N_min+1)), list(range(NJ))):
        N_, J_ = N[i], J[j]
        p.apply_async(func, args=(N_, J_, i, j, Phase,))

    p.close()
    p.join()
    P = [list(x) for x in list(Phase)]
    stime = nt()
    with open("PD_Data_"+stime+".txt", "w") as f:
        d = dict(N_min=N_min, N_max=N_max, J_min=J_min,
                 J_max=J_max, NJ=NJ, Phase=P)
        f.write(dumps(d))

    Draw(N, J, P, title=title,
         filename="PhaseDiag_"+stime)
    T_end = datetime.now()
    print("End Calculation at ", str(T_end))
    print("Total time:", str(T_end-T_start))


if __name__ == "__main__":
    # CalcGap()
    PhaseDiag(Run_4, title="Phase Diagram of N & J\n(Spin z in top and bottom layer)")
