import z2pack
import numpy as np
import os
from matplotlib import pyplot as plt

identity = np.identity(2, dtype=complex)
pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
pauli_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
pauli_z = np.array([[1, 0], [0, -1]], dtype=complex)
pauli_vector = list([pauli_x, pauli_y, pauli_z])

settings = {'num_lines': 50,
            'pos_tol': 1e-2,
            'gap_tol': 2e-2,
            'move_tol': 0.3,
            'iterator': range(50, 101, 2),
            'min_neighbour_dist': 1e-3,
            }


def h_QAH_p(M, B):
    r'''
    2*2 Hamiltonian for QAH;
    para A,M,B
    '''
    A = 1

    def hamiltonian(k):
        # 在z2pack中使用的是约化的k:  (kx,ky) in ([0,1],[0,1])
        # 在hamiltonian中需要将其恢复？
        kx, ky = k[0]*2*np.pi, k[1]*2*np.pi
        # kx, ky = k[0]*10, k[1]*10
        # 对于k·p model,由于它只在gamma点附近是准的，因此任意参数都是可以的
        return A*kx*pauli_x+A*ky*pauli_y+(M-B*(kx*kx+ky*ky))*pauli_z
    return hamiltonian


h0 = h_QAH_p(0.01, 1)
s0 = z2pack.hm.System(h0, dim=2, bands=1)
res1 = z2pack.surface.run(
    system=s0,
    # parameter of surface is moduled by 2pi
    surface=lambda k1, k2: [k2-0.5, k1/2-0.5],
    # surface=lambda k1, k2: [k2, k1/2],
    **settings
    # save_file="savefile.msgpack"
)
res2 = z2pack.surface.run(
    system=s0,
    # parameter of surface is moduled by 2pi
    surface=lambda k1, k2: [k2-0.5, k1/2],
    # surface=lambda k1, k2: [k2, k1/2],
    **settings
    # save_file="savefile.msgpack"
)
result = z2pack.surface.run(
    system=s0,
    # parameter of surface is moduled by 2pi
    surface=lambda k1, k2: [k2-0.5, k1-0.5],
    # surface=lambda k1, k2: [k2, k1/2],
    **settings
    # save_file="savefile.msgpack"
)
# print("Chern=", z2pack.invariant.chern(result))
print("Z2_1=", z2pack.invariant.z2(res1))
print("Z2_2=", z2pack.invariant.z2(res2))
print("Z2=", z2pack.invariant.z2(result))


fig, ax = plt.subplots(1,3)
z2pack.plot.wcc(res1, axis=ax[0])
z2pack.plot.wcc(res2, axis=ax[1])
z2pack.plot.wcc(result, axis=ax[2])

# z2pack.plot.wcc(res, axis=ax[1])
# ax.set_xlabel(r"$k_x$")
# ax.set_ylabel(r"$\bar{x}$")
# ax.set_title("QAH("+r"$\bf{k}\cdot\bf{p}$"+" model)")
# plt.show()
plt.savefig("QAH-Z2test.png")
