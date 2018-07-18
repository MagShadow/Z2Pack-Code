import numpy as np
import z2pack
import os
from matplotlib import pyplot as plt

# defining pauli matrices
identity = np.identity(2, dtype=complex)
pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
pauli_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
pauli_z = np.array([[1, 0], [0, -1]], dtype=complex)


def Hamilton(k, m, t1, t2, phi):
    # See "Topological Insulators", Shun-Qing Shen, Page 87
    kx, ky, _ = k
    k_a = 2 * np.pi / 3. * np.array([
        kx + ky,
        -2. * kx + ky,
        kx - 2. * ky
    ])
    k_b = 2 * np.pi * np.array([-kx + ky, -ky, kx])
    H = 2 * t2 * np.cos(phi) * sum(np.cos(k_b)) * identity
    H += t1 * sum(np.cos(k_a)) * pauli_x
    H += t1 * sum(np.sin(k_a)) * pauli_y
    H += m * pauli_z
    H -= 2 * t2 * np.sin(phi) * sum(np.sin(k_b)) * pauli_z
    return H


def get_chern(m, t1, t2, phi, **settings):
    # This function should return the Chern number for the given
    # parameters m, t1, t2, phi. The **settings should be passed
    # on to the z2pack.surface.run method.
    def h(k):
        return Hamilton(k, m, t1, t2, phi)

    # def honeycomb(k1, k2): return [(k1+k2)*np.sqrt(3)/2., (k1-k2)/2., 0]
    sys = z2pack.hm.System(h, bands=1)
    res = z2pack.surface.run(system=sys, surface=lambda k1, k2: [
                             k1, k2, 0], **settings)
    fig, ax = plt.subplots(1, 2)
    z2pack.plot.chern(res, axis=ax[0])
    z2pack.plot.wcc(res, axis=ax[1])
    plt.savefig("HaldaneModel.png")
    # print("Z2 = ", z2pack.invariant.z2(res))
    return z2pack.invariant.chern(res)


if __name__ == "__main__":
    # Task a)
    print(get_chern(0.5, 1, 1/3, 0.5*np.pi, num_lines=101,
                    min_neighbour_dist=1e-4))
    # print(get_chern(0.5, 1, 1/3, -0.5*np.pi))

    # Task b)

    # Task c) (optional - for advanced Python users)
