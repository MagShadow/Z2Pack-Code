# Compare the result in diffrent parameters, namely N, Delta & Constants

import numpy as np
from PhaseDiagram import PhaseDiag
from _utils import nt, convertSpin, settings, settings_strict, CONST_HJZ, CONST_HZL
from TI_Film import Hamiltonian
from TopoInvCalc import Calc, TopoOrder


def Run_0(_N, _J, _i, _j, Phase):
    '''
    This func calculate the system with uniform spin distribution in Z-direction.
    '''
    h = Hamiltonian(N=_N, J=_J, Delta=3, **CONST_HJZ)
    res = Calc(h, CalcZ2=True, LogOut=False, settings=settings)
    Phase[_i][_j] = TopoOrder(res)
    return


def Run_1(_N, _J, _i, _j, Phase):
    '''
    This func calculate the system with uniform spin distribution in Z-direction.
    '''
    h = Hamiltonian(N=_N, J=_J, Delta=3.333, **CONST_HJZ)
    res = Calc(h, CalcZ2=True, LogOut=False, settings=settings)
    Phase[_i][_j] = TopoOrder(res)
    return


now = nt()
File_0="PD_Delta_3_"+now
File_1="PD_Delta_3.333_"+now

PD_0 = PhaseDiag().run(Run_0, 6, 20, 15, 0, 0.02, 20, "N", "J")
PD_0.write(filename=File_0)
PD_0.draw(title=r"PhaseDiag of N & J, Delta=$3\AA$, Spin-z",
          xlabel="# of Layers", ylabel=r"J(\rm eV)", filename=File_0)

# File_1 = "PhaseDiag_"+nt()
PD_1 = PhaseDiag().run(Run_0, 6, 18, 13, 0, 0.02, 20, "N", "J")
PD_1.write(filename=File_1)
PD_1.draw(title=r"PhaseDiag of N & J, Delta=$3.33\AA$, Spin-z",
          xlabel="# of Layers", ylabel=r"J(\rm eV)", filename=File_1)
# PD_1 = PhaseDiag().run(Run_1, 6, 18, 13, 0, 0.02, 20, "N", "J")
