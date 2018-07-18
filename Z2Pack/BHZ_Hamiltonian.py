import z2pack
import numpy as np
import os
from matplotlib import pyplot as plt

identity = np.identity(2, dtype=complex)
pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
pauli_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
pauli_z = np.array([[1, 0], [0, -1]], dtype=complex)
pauli_vector = list([pauli_x, pauli_y, pauli_z])

settings = {'num_lines': 201,
            'pos_tol': 1e-2,
            'gap_tol': 2e-2,
            'move_tol': 0.3,
            'iterator': range(100, 201, 2),
            'min_neighbour_dist': 1e-3,
            }

# d=58 angstrom, A=-3.62, B=-18.0, C=-0.0180, D=-0.594, M=0.00922
# d=70 angstrom, A=-3.42, B=-16.9, C=-0.0263, D=0.514, M=-0.00686


def BHZ(M, B):
    # A, C, D = -3.62, -0.0180, -0.594
    A, C, D = -3.42, -0.0263, 0.514
    # A, C, D = 1, 0, 0
    #A, C, D = 1, 0, 1/2
    # k=(kx,ky,kz)
    # \epsilon(k)=C-D*k^2

    def h(k):
        kx, ky = k[0]*2*np.pi, k[1]*2*np.pi
        return (C-D*(kx*kx+ky*ky))*identity+A*kx*pauli_x+A*ky*pauli_y+(M-B*(kx*kx+ky*ky))*pauli_z
    # print(h([0, 0, 0]))

    def hamiltonian(k0):
        k = np.array(k0)
        mk = np.array([-x for x in k])
        # print(k0,mk)
        return np.column_stack((np.row_stack((h(k), np.zeros((2, 2)))), np.row_stack((np.zeros((2, 2)), h(mk).conjugate()))))

    return hamiltonian


# h0 = BHZ(0.00922, -18.0)
# h0 = BHZ(-0.00686, -16.9)
h0 = BHZ(-0.0686, -16.9)
# print(h0([0, 0, 0]))
s0 = z2pack.hm.System(h0, dim=2)
# sph = z2pack.shape.Sphere([0, 0, 0], 0.01)

# Only Consider Half


def squ(k1, k2): return np.array([k2-0.5, k1/2-0.5])
# def squ(k1, k2): return np.array([k2, k1/2])


result = z2pack.surface.run(
    system=s0,
    surface=squ,  # parameter of surface is moduled by 2pi
    **settings
    # save_file="savefile.msgpack"
)
print("Chern=", z2pack.invariant.chern(result))
print("Z2=", z2pack.invariant.z2(result))

fig, ax = plt.subplots()
# z2pack.plot.chern(result, axis=ax[0])
z2pack.plot.wcc(result, axis=ax)
# ax.set_xlim(0.8, 1.01)
# ax.set_xlim(-0.01, 0.2)
# ax[0].set_xlim(-1, 2)
# ax[0].set_ylim(-1, 2)
# ax[1].set_xlim(-1, 2)
# ax[1].set_ylim(-1, 2)
plt.savefig("BHZ.png")
