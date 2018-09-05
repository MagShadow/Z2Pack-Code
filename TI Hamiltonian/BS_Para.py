import numpy as np
# import matplotlib as mpl
# # mpl.use("Agg")
# from matplotlib import pyplot as plt
from PhaseDiagram import PhaseDiag
from _utils import nt, convertSpin, settings, settings_strict, CONST_HJZ, CONST_HZL
from TI_Film import Hamiltonian, plotLine
from TopoInvCalc import Calc, TopoOrder

if __name__ == '__main__':
    N, J = 25, 0.00
    h = Hamiltonian(N=N, J=J, Delta=1.9134, **CONST_HZL)
    res = Calc(h, CalcZ2=False, KScale=2, settings=settings_strict)
    res.plotChern(title="N=25 (5 QLs), J=0.00 eV, Spin-Z",
                  filename="HZL_25_0.00_1.9134_K2")
    # xRange, yRange, Nx, Ny = 0.05, 0.05, 50, 50
    # plotLine(h, 2*N-4, 2*N+4, xRange=xRange, Nx=Nx,
    #          title="Band Structure(8 bands), N=30, Delta=$1.9134\AA$, No Spin\nParameters from Hai-Zhou Lu's paper.", filename="N30_J0_D1.9134")
