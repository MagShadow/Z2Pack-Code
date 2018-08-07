# This program is aimed to search the Topological Phase
# in 2D phase diagram of N(hybrdization gap) & J(Spin Gap)

import numpy as np
import os
from scipy import linalg
from matplotlib import pyplot as plt

from TI_Film import Hamiltonian, Eig


def CalcGap():
    N_min, N_max, J = 6, 30, 0.00
    Gap = np.zeros([N_max+1], dtype=float)
    for N in range(N_min, N_max+1):
        S_ = np.array([[1, 0, 0]]*N)
        S = np.array([([s[0]*np.sin(s[1])*np.cos(s[2]), s[0]*np.sin(s[1])
                        * np.sin(s[2]), s[0]*np.cos(s[1])])for s in S_])
        h = Hamiltonian(N=N, J=J)
        eig = np.array([x.real for x in (linalg.eig(h(0, 0))[0])])
        eig.sort()
        Gap[N] = eig[2*N]-eig[2*N-1]

    plt.subplot()
    plt.plot(list(range(N_min, N_max+1)), Gap[N_min:N_max+1]*1000)
    plt.xlabel(r"$N$")
    plt.ylabel(r"$Gap(\rm meV)$")
    plt.title("Gap vs Thickness, No Spin")
    plt.savefig("Gap_vs_Thickness_No_Spin.png")


if __name__ == "__main__":
    CalcGap()
