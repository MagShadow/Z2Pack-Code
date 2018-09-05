import numpy as np

# 在服务器端没有DISPLAY资源，因此必须换成无交互的后端。默认为"TkAgg"。
import matplotlib as mpl
mpl.use("Agg")
from matplotlib import pyplot as plt

from itertools import product
from functools import partial
from scipy import linalg

from PhaseDiagram import PhaseDiag
from _utils import nt, convertSpin, settings, settings_strict, CONST_HJZ, CONST_HZL
from TI_Film import Hamiltonian, plotLine
from TopoInvCalc import Calc, TopoOrder
from runfunc import Run_SpinZ


def read(filename="PD_HZL_Delta_3_18-08-26-17-07-00.txt"):
    PD = PhaseDiag()
    PD.read(filename)
    X_L, X_H, X_N, X_Name = PD.info["X_L"], PD.info["X_H"], PD.info["X_N"], PD.info["X_Name"]
    Y_L, Y_H, Y_N, Y_Name = PD.info["Y_L"], PD.info["Y_H"], PD.info["Y_N"], PD.info["Y_Name"]
    X = np.linspace(X_L, X_H, num=X_N, endpoint=True)
    Y = np.linspace(Y_L, Y_H, num=Y_N, endpoint=True)
    with open("Res_HZL_O_18-08-26-17-07-00.txt", "w") as f:
        for i, j in product(list(range(X_N)), list(range(Y_N))):
            _X, _Y = X[i], Y[j]
            if PD.data[i][j] != 0:
                f.write("N=%d, J=%.3f; Res=%d\n" % (X[i], Y[j], PD.data[i][j]))
    return PD


def read2(filename="Poster_D3.189_SpinZ.txt"):
    PD = PhaseDiag()
    PD.read(filename)
    X_L, X_H, X_N, X_Name = PD.info["X_L"], PD.info["X_H"], PD.info["X_N"], PD.info["X_Name"]
    Y_L, Y_H, Y_N, Y_Name = PD.info["Y_L"], PD.info["Y_H"], PD.info["Y_N"], PD.info["Y_Name"]
    # Y_L, Y_H, Y_N = 0, 0.02, 11
    X = np.linspace(X_L, X_H, num=X_N, endpoint=True)
    Y = np.linspace(Y_L, Y_H, num=Y_N, endpoint=True)
    for i, j in product(list(range(X_N)), list(range(Y_N))):
        PD.data[i][j] = min(PD.data[i][j], 2)
    return PD


def read3(filename="Poster_PT_D3.189_HZL_SpinZ18-09-02-02-08-42.txt"):
    PD = PhaseDiag()
    PD.read(filename)
    X_L, X_H, X_N, X_Name = PD.info["X_L"], PD.info["X_H"], PD.info["X_N"], PD.info["X_Name"]
    Y_L, Y_H, Y_N, Y_Name = PD.info["Y_L"], PD.info["Y_H"], PD.info["Y_N"], PD.info["Y_Name"]
    # Y_L, Y_H, Y_N = 0, 0.02, 11
    X = np.linspace(X_L, X_H, num=X_N, endpoint=True)
    Y = np.linspace(Y_L, Y_H, num=Y_N, endpoint=True)
    # for i, j in product(list(range(X_N)), list(range(Y_N))):
    #     PD.data[i][j] = min(PD.data[i][j], 2)
    return [l[0] for l in PD.data]


def draw1():
    height = 0.04

    Filename = "Poster_PT"
    PT = read3()
    PT = np.array([x if (x > 0) and (x < height) else height for x in PT])
    X = np.linspace(3, 20, 18, endpoint=True, dtype=int)
    X1 = np.linspace(3, 19, 9, endpoint=True, dtype=int)
    X2 = np.linspace(4, 20, 9, endpoint=True, dtype=int)
    Y1 = [PT[i-3] for i in X1]
    Y2 = [PT[i-3] for i in X2]
    Y = np.array([height-x for x in PT])
    # print(X, Y)
    fig, ax = plt.subplots()
    plt.bar(X1, Y1, width=1, facecolor="green")
    plt.bar(X2, Y2, width=1, facecolor="red")
    plt.bar(X, Y, width=1, bottom=PT, facecolor="blue")
    plt.xlim(2.5, 20.5)
    plt.ylim(0, height)
    plt.title("PhaseDiag of N & J, Delta=$3.189\AA$")
    plt.xlabel("# of Layers")
    plt.ylabel("$J(\\rm eV)$")
    cmap = mpl.colors.ListedColormap(["r", "g", "b", "c"])
    norm = mpl.colors.BoundaryNorm(list(range(5)), cmap.N)
    ax2 = fig.add_axes([0.92, 0.1, 0.03, 0.8])
    mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm)
    # plt.show()
    plt.savefig("PT_D3.189.png")


def CalcGap():
    TEST_HJZ = {
        "A1": 2.26, "A2": 3.33, "C": -0.0083, "D1": 5.74, "D2": 30.4,
        "M": 0.28, "B1": 6.86, "B2": 44.5
    }
    TEST_HZL = {  # Interesting!!!!!
        "A1": 3.3, "A2": 4.1, "C": -0.0068, "D1": 1.2, "D2": -30.1,
        "M": 0.28, "B1": 1.5, "B2": -54.1
    }

    N_min, N_max, J = 5, 30, 0.00
    Gap = np.zeros([N_max-N_min+1], dtype=float)
    for N in range(N_min, N_max+1):
        # By Default it is Spin-Z
        h = Hamiltonian(N=N, J=J, Delta=1.9134, **TEST_HZL)
        eig = np.array([x.real for x in (linalg.eig(h(0, 0))[0])])
        eig.sort()
        Gap[N-N_min] = eig[2*N]-eig[2*N-1]

    plt.subplot()
    plt.plot(list(range(N_min, N_max+1)), Gap*1000)
    plt.xlabel(r"$N$")
    plt.ylabel(r"$Gap(\rm meV)$")
    plt.title("Gap vs Thickness, No Spin")
    plt.savefig("Gap_vs_Thickness_No_Spin_D1.9134.png")
    return


if __name__ == '__main__':
    # CalcGap()
    now = nt()
    File = "PD_SpinZ D1.9134_HJZ"+now
    # PD = read2("PD_HZL_Inv_N15_18-09-02-14-27-26.txt")
    # PD.draw(title="PhaseDiag of $\\theta$ & J, N=15, Delta=$3.189\AA$, Inversion",
    #         xlabel="$\\theta$", ylabel="$J(\\rm eV)$", filename=File_2)
    # print(PT, Y1, Y2)
    # PD.draw(title="PhaseDiag of N & J, Delta=$3.189\AA$, Spin-z",
    #         xlabel="# of Layers", ylabel="$J(\\rm eV)$", filename=Filename)

    func = partial(Run_SpinZ, Delta=1.9134, CONST=CONST_HJZ, settings=settings)
    # func2 = partial(Run_SpinZ, Delta=3.189, settings=settings)

    PD = PhaseDiag().run(func, 15, 30, 4, 0, 0.02, 3, "N", "J")
    # PD_2 = PhaseDiag().run(func2, 3, 20, 18, 0, 0.03, 16, "N", "J")
    # PD_2.write(filename=File_2)
    # PD_1 = PhaseDiag().run(func1, 3, 20, 18, 0, 0.03, 16, "N", "J")
    PD.write(filename=File)

    PD.draw(title="PhaseDiag of N & J, Delta=$1.9134\AA$, Spin-z\nParameters from Hai-Zhou Lu's paper.",
            xlabel="# of Layers", ylabel="$J(\\rm eV)$", filename=File)
    # PD_2.draw(title="PhaseDiag of N & J, Delta=$3.189\AA$, Spin-z\nParameters from Hai-Zhou Lu's paper.",
    #           xlabel="# of Layers", ylabel="$J(\\rm eV)$", filename=File_2)
