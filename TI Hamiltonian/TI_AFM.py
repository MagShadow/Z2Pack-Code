import os
import numpy as np
from scipy import linalg
from numbers import Integral
import matplotlib as mpl
# mpl.use("Agg")
from matplotlib import pyplot as plt

from TI_Film import Eig, plotLine
# from TI_Film import Hamiltonian as Ham_Flim
settings = {'num_lines': 81,
            'pos_tol':  5e-4,
            'gap_tol': 0.01,
            'move_tol': 0.4,
            'iterator': range(80, 121, 5),
            'min_neighbour_dist': 1e-5,
            }
identity = np.identity(2, dtype=complex)
pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
pauli_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
pauli_z = np.array([[1, 0], [0, -1]], dtype=complex)
pauli_vector = list([pauli_x, pauli_y, pauli_z])

# Some Constant
# use eV and Angstrom as units
m0, m1, hvf, t = 0.05, 16.45, 3.29, 0.005
hbar = 1.5
a_lat = 15
M_T, M_B = 0.1, 0.1


def Hamiltonian(M_T=M_T, M_B=M_B):
    def _m(kx, ky):
        return m0-m1*(kx*kx+ky*ky)

    def _h0(kx, ky):
        return hbar*hbar*(kx*kx+ky*ky)/m0

    def _ham(kx, ky):
        H = np.zeros([8, 8], dtype=complex)
        H[0:2, 0:2] = H[6:8, 6:8] = _h0(kx, ky)*identity
        H[2:4, 4:6] = H[4:6, 2:4] = _m(kx, ky)*identity
        H[2:4, 0:2] = H[6:8, 4:6] = t*identity
        H[0:2, 2:4] = H[4:6, 6:8] = np.conj(t)*identity
        H[2:4, 2:4] = hvf*(pauli_x*kx-pauli_y*ky)+M_T*pauli_z
        H[4:6, 4:6] = -hvf*(pauli_x*kx-pauli_y*ky)+M_B*pauli_z
        return H

    return _ham


def Ham_Small(M_T=M_T, M_B=M_B):
    def _m(kx, ky):
        return m0-m1*(kx*kx+ky*ky)

    def _h0(kx, ky):
        return hbar*hbar*(kx*kx+ky*ky)/m0

    def _ham(kx, ky):
        H = np.zeros([4, 4], dtype=complex)
        H[0:2, 0:2] = hvf*(pauli_x*kx-pauli_y*ky)+M_T*pauli_z
        H[2:4, 2:4] = -hvf*(pauli_x*kx-pauli_y*ky)+M_B*pauli_z
        H[0:2, 2:4] = H[2:4, 0:2] = _m(kx, ky)*identity
        # print(H)
        return H
    return _ham


if __name__ == "__main__":
    # h = Hamiltonian(M_T, M_B)
    # plotLine(h, 0, 8, xRange=2/a_lat, Nx=40)
    from TopoInvCalc import Calc
    # TRange, BRange, N_T, N_B = 0.1, 0.1, 10, 10
    # T = np.linspace(-TRange, TRange, num=N_T, endpoint=True)
    # B = np.linspace(-BRange, BRange, num=N_B, endpoint=True)
    t, b = 0.04, 0.06
    # print(t, b)
    h = Ham_Small(t, b)
    res = Calc(h, KScale=1, CalcZ2=False, LogOut=True, settings=settings)
    # print(res.Chern, res._Z2)
    res.plotChern()
    # plotLine(h, 0, 4, xRange=1/a_lat, Nx=40)
