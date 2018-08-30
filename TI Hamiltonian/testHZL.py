import numpy as np
from itertools import product
from functools import partial

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


if __name__ == '__main__':
    now = nt()
    File_1 = "PD_HZL_D3_Strict_"+now
    File_2 = "PD_HZL_D3.333_Strict_"+now
    func1 = partial(Run_SpinZ, Delta=3)
    func2 = partial(Run_SpinZ, Delta=3.333)

    PD_1 = PhaseDiag().run(func1, 3, 20, 18, 0, 0.03, 16, "N", "J")
    PD_1.write(filename=File_1)
    PD_2 = PhaseDiag().run(func2, 3, 20, 18, 0, 0.03, 16, "N", "J")
    PD_2.write(filename=File_2)

    PD_1.draw(title="PhaseDiag of N & J, Delta=$3\AA$, Spin-z\nParameters from Hai-Zhou Lu's paper.(Strict)",
              xlabel="# of Layers", ylabel="$J(\\rm eV)$", filename=File_1)
    PD_2.draw(title="PhaseDiag of N & J, Delta=$3.333\AA$, Spin-z\nParameters from Hai-Zhou Lu's paper.(Strict)",
              xlabel="# of Layers", ylabel="$J(\\rm eV)$", filename=File_2)
