import z2pack
import numpy as np
import os
from matplotlib import pyplot as plt

pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
pauli_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
pauli_z = np.array([[1, 0], [0, -1]], dtype=complex)
pauli_vector = list([pauli_x, pauli_y, pauli_z])

# d=70 angstrom, A=-3.42, B=-16.9, C=-0.0263, D=0.514, M=-0.00686


def BHZ(M, B):
    A, C, D = -3.42, -0.0263, 0.514
    #A, C, D = 1, 0, 1/2
    # k=(kx,ky,kz)
    # \epsilon(k)=C-D*k^2

    def h(k):
        return (C-D*(k[0]*k[0]+k[1]*k[1]))*np.identity(2)+A*k[0]*pauli_x+A*k[1]*pauli_y+(M-B*(k[0]*k[0]+k[1]*k[1]))
    # print(h([0,0,0]))
    def hamiltonian(k):
        return np.column_stack((np.row_stack((h(k), np.zeros((2, 2)))), np.row_stack((np.zeros((2, 2)), h(-k).conjugate()))))

    return hamiltonian


h0 = BHZ(-0.00686, -16.9)

