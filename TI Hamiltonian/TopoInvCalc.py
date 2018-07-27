import z2pack
import numpy as np
from matplotlib import pyplot as plt
import os

from TI_Film import Hamiltonian, Eig, plotBS
# Constant for Z2Pack Calculation
identity = np.identity(2, dtype=complex)
pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
pauli_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
pauli_z = np.array([[1, 0], [0, -1]], dtype=complex)
pauli_vector = list([pauli_x, pauli_y, pauli_z])

settings = {'num_lines': 31,
            'pos_tol':  1e-2,
            'gap_tol': 0.2,
            'move_tol': 0.05,
            'iterator': range(30, 201, 4),
            'min_neighbour_dist': 6e-5,
            }


def Calc(ham, surf=lambda k1, k2: [k1-0.5, k2-0.5], CalcZ2=True):
    KScale = 1

    def h0(k): return ham(k[0]/KScale, k[1]/KScale)
    s0 = z2pack.hm.System(h0, dim=2)
    result = z2pack.surface.run(
        system=s0,
        # parameter of surface is moduled by 2pi
        surface=surf,
        # surface=lambda k1, k2: [k2, k1/2],
        **settings
        # save_file="savefile.msgpack"
    )

    class Res(object):
        def __init__(self):
            self.Chern = z2pack.invariant.chern(result)
            if CalcZ2 == True:
                result_l = z2pack.surface.run(
                    system=s0, surface=lambda k1, k2: surf(k1/2, k2), **settings)
                result_r = z2pack.surface.run(
                    system=s0, surface=lambda k1, k2: surf(k1/2+0.5, k2), **settings)
                Z2_1, Z2_2, Z2_T = z2pack.invariant.z2(result_l), z2pack.invariant.z2(
                    result_r), z2pack.invariant.z2(result)
                self.Z2 = True if (Z2_1*Z2_2*(1-Z2_T) == 1) else False
                self._Z2 = [Z2_1, Z2_2, Z2_T]

        def plotChern(self, title="", filename=""):
            fig, ax = plt.subplots()
            z2pack.plot.wcc(result, axis=ax)
            ax.set_xlabel(r"$k_x$")
            ax.set_ylabel(r"$\bar{x}$")
            # ax.set_xlim(0.5, 1)
            ax.set_title(title+"\nChern="+str(self.Chern))
            # # plt.savefig("Chern TI film z spin.png")
            # ax.set_title("TI film, d=6nm, j=0.02")
            if filename!="":
                path=os.path.join("Pictures",filename+".png")
                plt.savefig(path)
            else:
                plt.show()
    return(Res())


if __name__ == "__main__":
    N, J = 20, 0.02
    S_ = np.zeros([N, 3])
    for i in range(N):
        S_[i, 0] = 1
        S_[i, 1] = 0
        S_[i, 2] = 0
    S = np.array([([s[0]*np.sin(s[1])*np.cos(s[2]), s[0]*np.sin(s[1])
                    * np.sin(s[2]), s[0]*np.cos(s[1])])for s in S_])
    h = Hamiltonian(N=N, J=J, S=S)

    res = Calc(h, CalcZ2=False)
    print(res.Chern)
    res.plotChern(title="Spin-z, J=0.02")
