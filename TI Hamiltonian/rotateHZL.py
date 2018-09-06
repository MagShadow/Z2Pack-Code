import numpy as np
from itertools import product
from functools import partial

from PhaseDiagram import PhaseDiag
from _utils import nt, convertSpin, settings, settings_strict, CONST_HJZ, CONST_HZL
from TI_Film import Hamiltonian, plotLine
from TopoInvCalc import Calc, TopoOrder
from runfunc import Run_SpinZ, Run_Mirror, Run_Inv, Run_M_I

if __name__ == '__main__':
    # now = nt()
    # File_1 = "PD_FS_SpinZ_N25_"+now
    # # File_2 = "PD_HZL_Inv_N15_"+now
    # # File_3 = "PD_HZL_M+I_N15_"+now

    # func1 = partial(Run_Mirror, N=25, Delta=1.9134, CONST=CONST_HZL)
    # # func2 = partial(Run_Inv, N=15, CONST=CONST_HZL)
    # # func3 = partial(Run_M_I, N=15, CONST=CONST_HZL)

    # PD_1 = PhaseDiag().run(func1, 0, np.pi, 13, 0, 0.025, 11, "Theta", "J")
    # PD_1.write(filename=File_1)
    # # PD_2 = PhaseDiag().run(func2, 0, np.pi, 13, 0, 0.02, 11, "Theta", "J")
    # # PD_2.write(filename=File_2)
    # # PD_3 = PhaseDiag().run(func3, 0, np.pi, 13, 0, 0.02, 11, "Theta", "J")
    # # PD_3.write(filename=File_3)
    File_1 = "PD_FS_SpinZ_N25_18-09-05-15-44-06"
    PD_1 = PhaseDiag()
    PD_1.read("PD_FS_SpinZ_N25_18-09-05-15-44-06.txt")
    PD_1.draw(title=r"PhaseDiag of $\theta$ & J, N=25 (5 QLs), Delta=$1.9134\AA$, Mirror",
              xlabel="r$\theta$", ylabel="$J(\\rm eV)$", filename=File_1)
    # PD_2.draw(title="PhaseDiag of Theta & J, N=15, Delta=$3.189\AA$, Inversion",
    #           xlabel="$\\theta$", ylabel="$J(\\rm eV)$", filename=File_2)
    # PD_3.draw(title="PhaseDiag of $\\theta$ & J, N=15, Delta=$3\AA$, Mirror+Inversion",
    #           xlabel="$\\theta$", ylabel="$J(\\rm eV)$", filename=File_3)
