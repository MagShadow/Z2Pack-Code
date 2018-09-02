# These file defines several "Run" funcs used in the
# calculation of Phase Diagram.

import numpy as np

from PhaseDiagram import PhaseDiag
from _utils import nt, convertSpin, settings, settings_strict, CONST_HJZ, CONST_HZL
from TI_Film import Hamiltonian
from TopoInvCalc import Calc, TopoOrder


def Run_SpinZ(_N, _J, _i, _j, Phase, Delta=3, CONST=CONST_HZL, settings=settings):
    '''
    This func calculate the system with uniform spin distribution in Z-direction.
    Use default constants from Hai-Zhou Lu's paper.
    '''
    print("Start Calc: N=%d, J=%.3f" % (int(_N), _J))
    h = Hamiltonian(N=int(_N), J=_J, Delta=Delta, **CONST)
    res = Calc(h, CalcZ2=True, LogOut=False, settings=settings)
    Phase[_i][_j] = TopoOrder(res)
    print("End Calc: N=%d, J=%.3f; Res=%d" % (int(_N), _J, TopoOrder(res)))

    return


def Run_Mirror(_T, _J, _i, _j, Phase, N=20, n=None, Delta=3, CONST=CONST_HZL):
    '''
    This func calculate the system with Mirror Symmetry,   which means that for the upper half and the lower half,the parallel components changes sign while the perpendicular term keep unchanged
    :para `n` defines how many layers (in upper half) contains the Spin term. In default, all of the upper half have the spin term.
    '''
    print("Start Calculation: Theta=%.3f , J=%.3f" % (_T, _J))
    _S = np.zeros([N, 3])
    layer = int(N/2) if n == None else min(int(n), int(N/2))
    for i in range(layer):
        _S[i, 0], _S[N-i-1, 0] = 1, 1
        _S[i, 1], _S[N-i-1, 1] = _T, _T
        _S[N-i-1, 2] = np.pi
    S = convertSpin(_S)

    h = Hamiltonian(N=N, J=_J, S=S, Delta=Delta, **CONST)
    res = Calc(h, CalcZ2=True, LogOut=False, settings=settings)
    Phase[_i][_j] = TopoOrder(res)

    print("End Calculation: Theta=%.3f pi , J=%.3f, Result: C=%.4f , Z2=%s" %
          (_T/np.pi, _J, res.Chern, str(res._Z2)))
    return


def Run_Inv(_T, _J, _i, _j, Phase, N=20, n=None, Delta=3, CONST=CONST_HZL):
    '''
    This func calculate the system with Mirror Symmetry,   which means that for the upper half and the lower half,the parallel components changes sign while the perpendicular term keep unchanged
    :para `n` defines how many layers (in upper half) contains the Spin term. In default, all of the upper half have the spin term.
    '''
    print("Start Calculation: Theta=%.3f , J=%.3f" % (_T, _J))
    _S = np.zeros([N, 3])
    layer = int(N/2) if n == None else min(int(n), int(N/2))
    for i in range(layer):
        _S[i, 0], _S[N-i-1, 0] = 1, 1
        _S[i, 1], _S[N-i-1, 1] = _T, _T
    S = convertSpin(_S)

    h = Hamiltonian(N=N, J=_J, S=S, Delta=Delta, **CONST)
    res = Calc(h, CalcZ2=True, LogOut=False, settings=settings)
    Phase[_i][_j] = TopoOrder(res)

    print("End Calculation: Theta=%.3f pi , J=%.3f, Result: C=%.4f , Z2=%s" %
          (_T/np.pi, _J, res.Chern, str(res._Z2)))
    return


def Run_M_I(_T, _J, _i, _j, Phase, N=20, n=None, Delta=3, CONST=CONST_HZL):
    '''
     This func calculate the system with Mirror + Inversion Symmetry, which means that for the upper half and the lower half,the perpendicular components changes sign while the parallel term keep unchanged.

    :para `n` defines how many layers (in upper half) contains the Spin term. In default, all of the upper half have the spin term.
    '''
    print("Start Calculation: Theta=%.3f , J=%.3f" % (_T, _J))
    _S = np.zeros([N, 3])
    layer = int(N/2) if n == None else min(int(n), int(N/2))
    for i in range(layer):
        _S[i, 0], _S[N-i-1, 0] = 1, 1
        _S[i, 1], _S[N-i-1, 1] = _T, np.pi-_T
    S = convertSpin(_S)

    h = Hamiltonian(N=N, J=_J, S=S, Delta=Delta, **CONST)
    res = Calc(h, CalcZ2=True, LogOut=False, settings=settings)
    Phase[_i][_j] = TopoOrder(res)

    print("End Calculation: Theta=%.3f pi , J=%.3f, Result: C=%.4f , Z2=%s" %
          (_T/np.pi, _J, res.Chern, str(res._Z2)))
    return
