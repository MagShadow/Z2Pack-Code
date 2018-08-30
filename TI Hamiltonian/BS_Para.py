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
    h = Hamiltonian(N, Delta=3.3333, **CONST_HZL)
    xRange, yRange, Nx, Ny = 0.05, 0.05, 50, 50
    plotLine(h, 2*N-4, 2*N+4, xRange=xRange, Nx=Nx,
             title="Band Structure(12 bands), Delta=$3.3333\AA$, No Spin\nParameters from Hai-Zhou Lu's paper.", filename="Test")
