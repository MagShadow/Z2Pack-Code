import numpy as np
# import matplotlib as mpl
# # mpl.use("Agg")
# from matplotlib import pyplot as plt
from PhaseDiagram import PhaseDiag
from _utils import nt, convertSpin, settings, settings_strict, CONST_HJZ, CONST_HZL
from TI_Film import Hamiltonian, plotLine
from TopoInvCalc import Calc, TopoOrder

if __name__ == '__main__':
    N, J = 15, 0.004
    h = Hamiltonian(N=N, J=J, Delta=3.189, **CONST_HZL)
    res = Calc(h, CalcZ2=False, KScale=5, settings=settings_strict)
    # print(res.result.wcc)
    res.plotChern(title="HZL N=15, J=0.004, Delta=$3.189\AA$ Spin-Z",
                  filename="HJZ_15_0.004_3.189_K5")
    # xRange, yRange, Nx, Ny = 0.05, 0.05, 50, 50
    # plotLine(h, 2*N-4, 2*N+4, xRange=xRange, Nx=Nx,
    #          title="Band Structure(8 bands), N=15, Delta=$3.189\AA$, No Spin\nParameters from Hai-Zhou Lu's paper.", filename="N15_J0_D3.189.png")
