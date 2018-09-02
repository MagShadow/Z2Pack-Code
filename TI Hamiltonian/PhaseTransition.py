import numpy as np
from itertools import product
from functools import partial

from PhaseDiagram import PhaseDiag
from _utils import nt, convertSpin, settings, settings_strict, CONST_HJZ, CONST_HZL
from TI_Film import Hamiltonian, plotLine
from TopoInvCalc import Calc, TopoOrder


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
    now = nt()
    File_1 = "PT_D3.189_HZL_SpinZ"+now
    func1 = partial(PT_SpinZ, _J_Min=0, _J_tol=1e-5,
                    Delta=3.189, settings=settings)

    PD_1 = PhaseDiag().run(func1, 3, 20, 18, 0.1, 0.1, 1, "N", "J")
    PD_1.write(filename=File_1)
