import z2pack
import numpy as np
import os
import matplotlib.pyplot as plt


pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
pauli_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
pauli_z = np.array([[1, 0], [0, -1]], dtype=complex)
pauli_vector = list([pauli_x, pauli_y, pauli_z])


def hamiltonian0(k):
    """simple 2-band hamiltonian k.sigma with a Weyl point at k=0"""
    res = np.zeros((2, 2), dtype=complex)
    for kval, p_mat in zip(k, pauli_vector):
        res += kval*p_mat
    return res


def hamiltonian1(k):
    """simple 2-band hamiltonian k.sigma with k_y -> -k_y"""
    k[2] = -k[2]
    res = np.zeros((2, 2), dtype=complex)
    for kval, p_mat in zip(k, pauli_vector):
        res += kval*p_mat
    return res


System0 = z2pack.hm.System(hamiltonian0, bands=1)
System1 = z2pack.hm.System(hamiltonian1)
sph = z2pack.shape.Sphere([0, 0, 0], 0.01)
res0 = z2pack.surface.run(system=System0, surface=sph, num_lines=21)
res1 = z2pack.surface.run(system=System1, surface=sph)

fig, ax = plt.subplots(1, 2)

# plot styling
fs = 15
ax[0].set_xlabel(r'$\theta$', fontsize=fs)
ax[0].set_ylabel(r'$\bar{x}$', rotation='horizontal', fontsize=fs)
ax[1].set_xlabel(r'$\theta$', fontsize=fs)
ax[0].set_xticks([0, 1])
ax[1].set_xticks([0, 1])
ax[0].set_xticklabels([r'$0$', r'$\pi$'])
ax[1].set_xticklabels([r'$0$', r'$\pi$'])
ax[0].set_title(r'$\vec{k}.\vec{\sigma}$', fontsize=fs)
ax[1].set_title(r'$(k_x, -k_y, k_z).\vec{\sigma}$', fontsize=fs)

# plotting the evolution of polarization
z2pack.plot.chern(res0, axis=ax[0])
# z2pack.plot.wcc(res0, axis=ax[1])
z2pack.plot.chern(res1, axis=ax[1])
plt.show()

# plt.savefig('plot.pdf', bbox_inches='tight')


# Isotropic Weyl Point
# def hamiltonian(k):
#     kx, ky, kz = k
#     return np.array([[kz, kx-1j*ky], [kx+1j*ky]], dtype=complex)


# system = z2pack.hm.System(hamiltonian)

# result = z2pack.surface.run(
#     system=system,
#     surface=lambda t1, t2: np.array([t1, t2, 0]),
#     # num_lines=101
#     # save_file=os.path.join("Topological Invariant","z2pack",'test_savefile.json')
# )

# print(z2pack.invariant.chern(result))
# print(z2pack.invariant.z2(result))
