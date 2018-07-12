import z2pack
import numpy as np
import os
from matplotlib import pyplot as plt

pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
pauli_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
pauli_z = np.array([[1, 0], [0, -1]], dtype=complex)
pauli_vector = list([pauli_x, pauli_y, pauli_z])


def h_QAH_p(M, B):
    r'''
    2*2 Hamiltonian for QAH;
    para A,M,B
    '''
    A = 1

    def hamiltonian(k):
        return A*k[0]*pauli_x+A*k[1]*pauli_y+(M-B*(k[0]*k[0]+k[1]*k[1]))*pauli_z
    return hamiltonian


h0 = h_QAH_p(1, 1)
s0 = z2pack.hm.System(h0, dim=2,bands=[1])
result = z2pack.surface.run(
    system=s0,
    surface=lambda k1, k2: [k1, k2],  # parameter of surface is moduled by 2pi
    num_lines=11
    # save_file="savefile.msgpack"
)

print(z2pack.invariant.chern(result))
fig, ax = plt.subplots(1, 2)
z2pack.plot.chern(result, axis=ax[0])
z2pack.plot.wcc(result, axis=ax[1])
plt.show()