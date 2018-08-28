import numpy as np
# import matplotlib as mpl
# # mpl.use("Agg")
# from matplotlib import pyplot as plt
from PhaseDiagram import PhaseDiag
from _utils import nt, convertSpin, settings, settings_strict, CONST_HJZ, CONST_HZL
from TI_Film import Hamiltonian, plotLine
from TopoInvCalc import Calc, TopoOrder

if __name__ == '__main__':
    N = 20
    h = Hamiltonian(N, **CONST_HJZ)
    xRange, yRange, Nx, Ny = 0.05, 0.05, 50, 50
    plotLine(h, 2*N-2, 2*N+2, xRange=xRange, Nx=Nx,
             title="Band Structure(4 bands), Delta=$3\AA$, No Spin\nParameters from Hai-Jun Zhang's paper.", filename="Test")
