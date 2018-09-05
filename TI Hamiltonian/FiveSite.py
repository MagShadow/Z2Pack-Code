import numpy as np

from itertools import product
from functools import partial
from scipy import linalg

from PhaseDiagram import PhaseDiag
from _utils import nt, convertSpin, settings, settings_strict, CONST_HJZ, CONST_HZL
from TI_Film import Hamiltonian, plotLine
from TopoInvCalc import Calc, TopoOrder
from runfunc import Run_SpinZ

if __name__ == '__main__':
    now = nt()
    File = "PD_D1.9134_HZL"+now
    func = partial(Run_SpinZ, Delta=1.9134, CONST=CONST_HZL, settings=settings)
    PD = PhaseDiag().run(func, 10, 35, 6, 0, 0.03, 7, "N", "J")
    PD.write(filename=File)
    PD.draw(title="PhaseDiag of N & J, Delta=$1.9134\AA$, Spin-z\nParameters from Hai-Zhou Lu's paper.",
            xlabel="# of Layers", ylabel="$J(\\rm eV)$", filename=File)
