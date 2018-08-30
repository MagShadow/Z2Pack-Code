# Compare the result in diffrent parameters, namely N, Delta & Constants

import numpy as np
from PhaseDiagram import PhaseDiag
from _utils import nt, convertSpin, settings, settings_strict, CONST_HJZ, CONST_HZL
from TI_Film import Hamiltonian
from TopoInvCalc import Calc, TopoOrder

from sys import stderr


def Run_0(_N, _J, _i, _j, Phase):
    '''
    This func calculate the system with uniform spin distribution in Z-direction.
    '''
    print("Start Calc: N=%d, J=%.3f" % (int(_N), _J))
    h = Hamiltonian(N=int(_N), J=_J, Delta=3, **CONST_HJZ)
    res = Calc(h, CalcZ2=True, LogOut=False, settings=settings)
    Phase[_i][_j] = TopoOrder(res)
    print("End Calc: N=%d, J=%.3f; Res=%d" % (int(_N), _J, TopoOrder(res)))

    return


def Run_1(_N, _J, _i, _j, Phase):
    '''
    This func calculate the system with uniform spin distribution in Z-direction.
    '''
    print("Start Calc: N=%d, J=%.3f" % (int(_N), _J))
    h = Hamiltonian(N=int(_N), J=_J, Delta=3.333, **CONST_HJZ)
    res = Calc(h, CalcZ2=True, LogOut=False, settings=settings)
    Phase[_i][_j] = TopoOrder(res)
    print("End Calc: N=%d, J=%.3f; Res=%d" % (int(_N), _J, TopoOrder(res)))
    return


def Run_HZL(_N, _J, _i, _j, Phase):
    '''
    This func calculate the system with uniform spin distribution in Z-direction.
    Use the Parameters from the Hai-Zhou Lu's paper.
    '''
    print("Start Calc: N=%d, J=%.3f" % (int(_N), _J))
    h = Hamiltonian(N=int(_N), J=_J, Delta=3, **CONST_HZL)
    res = Calc(h, CalcZ2=True, LogOut=False, settings=settings)
    Phase[_i][_j] = TopoOrder(res)

    print("End Calc: N=%d, J=%.3f; Res=%d" % (int(_N), _J, TopoOrder(res)))
    if TopoOrder(res) != 0:
        stderr.write("!!!!!!!!\nN=%d, J=%.3f; Res=%d\n!!!!!!!!!!\n" %
                     (int(_N), _J, TopoOrder(res)))
    return


if __name__ == '__main__':
    now = nt()
    # File_0 = "PD_Delta_3_"+now
    File_1 = "PD_HZL_Origin_"+now

    # PD_0 = PhaseDiag().run(Run_0, 6, 20, 15, 0, 0.02, 21, "N", "J")
    # PD_0.write(filename=File_0)
    # PD_0.draw(title=r"PhaseDiag of N & J, Delta=$3\AA$, Spin-z",
    #           xlabel="# of Layers", ylabel=r"J(\rm eV)", filename=File_0)
    # File_0 = "PD_Delta_3_18-08-26-17-07-00.txt"
    # File_1 = "PD_Delta_3.33_18-08-26-17-07-00.txt"
    # File_2 = "PD_HZL_Delta_3_18-08-26-17-07-00.txt"

    # PD_0 = PhaseDiag()
    # PD_0.read(File_0)
    # PD_0.draw(title=r"PhaseDiag of N & J, Delta=$3\AA$, Spin-z",
    #           xlabel="# of Layers", ylabel=r"J(\rm eV)", filename=File_0[:-4])
    # PD_1 = PhaseDiag().run(Run_1, 6, 18, 13, 0, 0.02, 21, "N", "J")
    # PD_1.write(filename=File_1)
    PD_2 = PhaseDiag().run(Run_HZL, 6, 20, 5, 0, 0.5, 6, "N", "J")
    PD_2.write(filename=File_1)

    # PD_1.draw(title=r"PhaseDiag of N & J, Delta=$3.33\AA$, Spin-z",
    #           xlabel="# of Layers", ylabel=r"J(\rm eV)", filename=File_1)
    PD_2.draw(title=r"PhaseDiag of N & J, Delta=$3\AA$, Spin-z\nParameters from Hai-Zhou Lu's paper.",
              xlabel="# of Layers", ylabel=r"$J(\rm eV)$", filename=File_1)
