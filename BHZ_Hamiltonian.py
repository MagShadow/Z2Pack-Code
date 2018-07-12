import z2pack
import numpy as np
import os
from matplotlib import pyplot as plt

pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
pauli_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
pauli_z = np.array([[1, 0], [0, -1]], dtype=complex)
pauli_vector = list([pauli_x, pauli_y, pauli_z])

# d=58 angstrom, A=-3.62, B=-18.0, C=-0.0180, D=-0.594, M=0.00922
# d=70 angstrom, A=-3.42, B=-16.9, C=-0.0263, D=0.514, M=-0.00686


def BHZ(M, B):
    # A, C, D = -3.62, -0.0180, -0.594
    A, C, D = -3.42, -0.0263, 0.514
    #A, C, D = 1, 0, 1/2
    # k=(kx,ky,kz)
    # \epsilon(k)=C-D*k^2

    def h(k):
        kx, ky = k[0], k[1]
        return (C-D*(kx*kx+ky*ky))*np.identity(2)+A*kx*pauli_x+A*ky*pauli_y+(M-B*(kx*kx+ky*ky))*pauli_z
    # print(h([0, 0, 0]))

    def hamiltonian(k0):
        k = np.array(k0)
        mk = np.array([-x for x in k])
        # print(k0,mk)
        return np.column_stack((np.row_stack((h(k), np.zeros((2, 2)))), np.row_stack((np.zeros((2, 2)), h(mk).conjugate()))))

    return hamiltonian


# h0 = BHZ(0.00922, -18.0)
h0 = BHZ(-0.00686, -16.9)
# print(h0([0, 1, 0]))
s0 = z2pack.hm.System(h0, dim=2, bands=[0,1])
# sph = z2pack.shape.Sphere([0, 0, 0], 0.01)


def squ(l1, l2): return np.array([l1, l2])


result = z2pack.surface.run(
    system=s0,
    surface=squ,  # parameter of surface is moduled by 2pi
    num_lines=101,
    # save_file="savefile.msgpack"
)
print(z2pack.invariant.chern(result))
print(z2pack.invariant.z2(result))

fig, ax = plt.subplots(1, 2)
z2pack.plot.chern(result, axis=ax[0])
z2pack.plot.wcc(result, axis=ax[1])
plt.show()
