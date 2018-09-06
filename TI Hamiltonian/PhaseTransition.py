import numpy as np
from itertools import product
from functools import partial
# 在服务器端没有DISPLAY资源，因此必须换成无交互的后端。默认为"TkAgg"。
import matplotlib as mpl
mpl.use("Agg")
from matplotlib import pyplot as plt

from PhaseDiagram import PhaseDiag
from _utils import nt, convertSpin, settings, settings_strict, CONST_HJZ, CONST_HZL
from TI_Film import Hamiltonian, plotLine
from TopoInvCalc import Calc, TopoOrder


def read_point(filename="Poster_PT_D3.189_HZL_SpinZ18-09-02-02-08-42.txt"):
    PD = PhaseDiag()
    PD.read(filename)
    return [l[0] for l in PD.data]


def draw(PT, height=0.04, filename="Poster_PT"):

    PT = np.array([x if (x > 0) and (x < height) else height for x in PT])
    X = np.linspace(2, 7, 6, endpoint=True, dtype=int)
    X1 = np.linspace(3, 7, 3, endpoint=True, dtype=int)
    X2 = np.linspace(2, 6, 3, endpoint=True, dtype=int)
    Y1 = [PT[i-2] for i in X1]
    Y2 = [PT[i-2] for i in X2]
    Y = np.array([height-x for x in PT])
    # print(X, Y)
    fig, ax = plt.subplots()
    plt.bar(X1, Y1, width=1, facecolor="green")
    plt.bar(X2, Y2, width=1, facecolor="red")
    plt.bar(X, Y, width=1, bottom=PT, facecolor="blue")
    plt.xlim(1.5, 7.5)
    plt.ylim(0, height)
    plt.title("PhaseDiag of Thickness & J, Delta=$1.9134\AA$")
    plt.xlabel("Thickness (QL)")
    plt.ylabel("$J(\\rm eV)$")
    cmap = mpl.colors.ListedColormap(["r", "g", "b", "c"])
    norm = mpl.colors.BoundaryNorm(list(range(5)), cmap.N)
    ax2 = fig.add_axes([0.92, 0.1, 0.03, 0.8])
    mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm)
    # plt.show()
    plt.savefig("PT_D1.9134.png")


def PT_SpinZ(_N, _J_Max, _i, _j, Phase, _J_Min=0,  _J_tol=1e-4, Delta=3.189, CONST=CONST_HZL, settings=settings):
    assert _J_tol > 0, "_J_tol should >0!"

    print("Start Calc: N=%d" % (int(_N)))
    L, R = _J_Min, _J_Max

    print("===>Calculate N=%d, J=%.3f" % (int(_N), L))
    h = Hamiltonian(N=int(_N), J=L, Delta=Delta, **CONST)
    res = Calc(h, CalcZ2=True, LogOut=False, settings=settings)
    T_L = TopoOrder(res)

    print("Calculate N=%d, J=%.3f" % (int(_N), R))
    h = Hamiltonian(N=int(_N), J=R, Delta=Delta, **CONST)
    res = Calc(h, CalcZ2=True, LogOut=False, settings=settings)
    T_R = TopoOrder(res)

    if T_L == T_R:
        Phase[_i][_j] = -1
        print("===>End Calc: N=%d. No phase transition, Res=%d" %
              (int(_N), TopoOrder(res)))
    else:
        while (R-L) > _J_tol:
            M = (R+L)/2
            print("Calculate N=%d, J=%.3f" % (int(_N), M))
            h = Hamiltonian(N=int(_N), J=M, Delta=Delta, **CONST)
            res = Calc(h, CalcZ2=True, LogOut=False, settings=settings)
            T_M = TopoOrder(res)
            if T_M == T_R:
                R = M
            elif T_M == T_L:
                L = M
            else:
                Phase[_i][_j] = R
                raise Exception("Middle result different with either!")
        Phase[_i][_j] = R
        print("===>End Calc: N=%d, phase transition at %.3f" %
              (int(_N), R))
    return


if __name__ == "__main__":
    # now = nt()
    # File_1 = "PT_D1.9134_HZL_SpinZ"+now
    # func1 = partial(PT_SpinZ, _J_Min=0, _J_tol=1e-5,
    #                 Delta=1.9134, settings=settings)

    # PD_1 = PhaseDiag().run(func1, 10, 35, 6, 0.1, 0.1, 1, "N", "J")
    # PD_1.write(filename=File_1)
    data = read_point("PT_D1.9134_HZL_SpinZ18-09-05-08-46-02.txt")

    draw(data, height=0.05, filename="PT_D1.9134")
