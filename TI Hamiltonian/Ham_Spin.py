# Bi2Se3 film + Spin Term

import numpy as np
from scipy import linalg
from matplotlib import pyplot as plt
import z2pack
import os
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

# Constants related to real material
# unit of length: angstrom
# unit of energy: eV
A1, A2, C, D1, D2 = 2.26, 3.33, -0.0083, 5.74, 30.4
M, B1, B2 = 0.28, 6.86, 44.5
Delta = 3
N = 20  # total layers in z direction, 0.3nm
J = 0.02  # an arbitary value for spin coupling energy
# J = 0.01

# S_[i] discribe the spin distribution
# first component is the intensity, from 0 to 1
# second and three term is the theta and phi, describe the direction of the spin
# S[i] is the spin vector

S_ = np.zeros([N, 3])
for i in range(N):
    S_[i, 0] = 1
    S_[i, 1] = np.pi/2
    S_[i, 2] = np.pi/2


S = np.array([([s[0]*np.sin(s[1])*np.cos(s[2]), s[0]*np.sin(s[1])
                * np.sin(s[2]), s[0]*np.cos(s[1])])for s in S_])
# print(S)

# U[i] discribe the potential distribution on each site
U = np.zeros([N])

# Prerequisite matrices
Gzz = np.array(np.diag([D1-B1, D1+B1, D1-B1, D1+B1]), dtype=complex)
Gz = np.zeros([4, 4], dtype=complex)
Gz[0, 1] = Gz[1, 0] = 1.0j*A1
Gz[3, 2] = Gz[2, 3] = -1.0j*A1
t_p = -1*Gzz/(Delta*Delta)-Gz/(2*Delta)
t_m = -1*Gzz/(Delta*Delta)+Gz/(2*Delta)
# sig_xyz = pauli_xyz tensor identity
sig_x, sig_y, sig_z = [np.zeros([4, 4], dtype=complex) for i in range(3)]
sig_x[0, 2], sig_x[1, 3], sig_x[2, 0], sig_x[3, 1] = 1, 1, 1, 1
sig_y[0, 2], sig_y[1, 3], sig_y[2, 0], sig_y[3, 1] = -1j, -1j, 1j, 1j
sig_z[0, 0], sig_z[1, 1], sig_z[2, 2], sig_z[3, 3] = 1, 1, -1, -1
sig = np.array([sig_x, sig_y, sig_z])

# print(sum([S[2, j]*sig[j] for j in range(3)]))
# print(t_p.T.conjugate()-t_m)
# print(t_m)
# print(S)
# Define a discretized Hamitonian


def Hamiltonian(kx, ky):
    def E_(kx, ky):
        return C+D2*(kx*kx+ky*ky)

    def M_(kx, ky):
        return M-B2*(kx*kx+ky*ky)

    h0 = np.zeros([4, 4], dtype=complex)
    h0[0, 0] = h0[2, 2] = M_(kx, ky)
    h0[1, 1] = h0[3, 3] = -1*M_(kx, ky)
    h0[0, 3] = h0[1, 2] = A2*(kx-1.0j*ky)
    h0[3, 0] = h0[2, 1] = A2*(kx+1.0j*ky)
    zero = np.zeros([4, 4], dtype=complex)
    # print(h0)
    # hii = E_(kx, ky)*np.identity(4, dtype=complex)+2 * \
    #     Gzz/(Delta*Delta)+h0
    h_ = np.array([((E_(kx, ky)+U[i])*np.identity(4, dtype=complex)+2 *
                    Gzz/(Delta*Delta)+h0-J*sum([S[i, j]*sig[j] for j in range(3)])) for i in range(N)])
    # print(hii)

    def row(i):
        return np.column_stack([(h_[i] if j == i else (t_m if j+1 == i else(t_p if j-1 == i else zero)))for j in range(N)])
    return np.row_stack([row(i) for i in range(N)])


# use Z2pack to calculate the Chern number
identity = np.identity(2, dtype=complex)
pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
pauli_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
pauli_z = np.array([[1, 0], [0, -1]], dtype=complex)
pauli_vector = list([pauli_x, pauli_y, pauli_z])

settings = {'num_lines': 101,
            'pos_tol': 5e-3,
            'gap_tol': 0.2,
            'move_tol': 0.05,
            'iterator': range(60, 201, 4),
            'min_neighbour_dist': 6e-5,
            }


# def h0(k): return Hamiltonian(k[0]/5, k[1]/5)


# s0 = z2pack.hm.System(h0, dim=2)
# result = z2pack.surface.run(
#     system=s0,
#     # parameter of surface is moduled by 2pi
#     surface=lambda k1, k2: [k1-0.5, k2-0.5],
#     # surface=lambda k1, k2: [k2, k1/2],
#     **settings
#     # save_file="savefile.msgpack"
# )
# print("Chern=", z2pack.invariant.chern(result))
# res1 = z2pack.surface.run(
#     system=s0,
#     # parameter of surface is moduled by 2pi
#     surface=lambda k1, k2: [k2-0.5, k1/2-0.5],
#     # surface=lambda k1, k2: [k2, k1/2],
#     **settings,
#     # save_file="savefile1.msgpack"
# )
# Z2_1= z2pack.invariant.z2(res1)
# print("Z2_1=",Z2_1)
# res2 = z2pack.surface.run(
#     system=s0,
#     # parameter of surface is moduled by 2pi
#     surface=lambda k1, k2: [k2-0.5, k1/2],
#     # surface=lambda k1, k2: [k2, k1/2],
#     **settings,
#     # save_file="savefile2.msgpack"
# )
# Z2_2 = z2pack.invariant.z2(res2)
# print("Z2_2=", Z2_2)

# result = z2pack.surface.run(
#     system=s0,
#     # parameter of surface is moduled by 2pi
#     surface=lambda k1, k2: [k2-0.5, k1-0.5],
#     # surface=lambda k1, k2: [k2, k1/2],
#     **settings,
#     save_file="savefile.msgpack"
# )
# Z2_1, Z2_2, Z2 = z2pack.invariant.z2(res1), z2pack.invariant.z2(
#     res2), z2pack.invariant.z2(result)
# print("Z2=", Z2_1, Z2_2, Z2)
# print("Chern=", z2pack.invariant.z2(result))

# fig, ax = plt.subplots(1, 3)
# z2pack.plot.wcc(res1, axis=ax[0])
# z2pack.plot.wcc(res2, axis=ax[1])
# z2pack.plot.wcc(result, axis=ax[2])
# # plt.savefig("BHZ-Z2test.png")
# ax[0].set_xlim(0.9, 1)
# ax[1].set_xlim(0, 0.1)
# ax[2].set_xlim(0.45, 0.55)

# ax[0].set_title("Z2_1="+str(Z2_1))
# ax[1].set_title("Z2_2="+str(Z2_2))
# ax[2].set_title("Z2_Total="+str(Z2))

# fig, ax = plt.subplots()
# z2pack.plot.wcc(res2, axis=ax)
# ax.set_xlabel(r"$k_x$")
# ax.set_ylabel(r"$\bar{x}$")
# # ax.set_xlim(0.5, 1)
# ax.set_title("TI film y spin, WCC2, j=0.06")
# # # plt.savefig("Chern TI film z spin.png")
# # ax.set_title("TI film, d=6nm, j=0.02")
# plt.savefig("WCC2 TI film y spin j=0.06.png")
# plt.show()


# # print(Hamiltonian(0.03, 0))
# # Simulation parameter, unit = 1/anstrom
# from -xRange to xRange, etc
xRange, yRange, Nx, Ny = 0.05, 0.05, 20, 20
dkx, dky = 2*xRange/Nx, 2*yRange/Ny
bs = np.zeros([Nx+1, Ny+1, 4*N], dtype=float)

# # do simulation on ky=0 line
# ky = 0
# j0 = j = int(Ny/2)
# for i in range(Nx+1):
#     # kx, ky = -xRange+dkx*i, -yRange+dky*j
#     kx = -xRange+dkx*i
#     # print(kx, end=", ")
#     temp = np.array([x.real for x in (linalg.eig(Hamiltonian(kx, ky))[0])])
#     temp.sort()
#     # print(temp)
#     bs[i, j] = temp
#     # print(temp)

# do simulation on kx-ky plane
for i in range(Nx+1):
    for j in range(Ny+1):
        kx, ky = -xRange+dkx*i, -yRange+dky*j
        temp = np.array([x.real for x in (linalg.eig(Hamiltonian(kx, ky))[0])])
        temp.sort()
        bs[i, j] = temp

Eig = np.array([([([bs[i, j, n] for j in range(Ny+1)]) for i in range(Nx+1)])
                for n in range(4*N)])

# print(Eig)
# print(Eig[2*N-1])
#####################################################
#3D plot
X = np.linspace(-xRange, xRange, Nx+1, endpoint=True)
Y = np.linspace(-yRange, yRange, Ny+1, endpoint=True)
X, Y = np.meshgrid(X, Y)
# print("X=", X)
# print("Y=", X)

fig = plt.figure()
ax = fig.gca(projection='3d')
for Z in Eig[2*N-2:2*N+2]:
    # print(Z)
    surf = ax.plot_surface(X, Y, Z,  cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
plt.show()
#####################################################
# j0 = int(Ny/2)
# X = np.linspace(-xRange, xRange, Nx+1, endpoint=True)
# # Eig = [([bs[j, j0, i] for j in range(Nx+1)])for i in range(4*N)]


# plt.subplot(1, 1, 1)
# for t in Eig[2*N-1:2*N]:
#     print(t)
#     y=np.array([t[i][j0] for i in range(Nx+1)])
#     plt.plot(X, y)
# plt.ylim(0.15, 0.35)
# plt.xlim(-0.05, 0.05)


# # plt.ylim(0.1, 0.5)
# # plt.xlim(-0.05, 0.05)
# plt.xlabel(r"$k_x(\rm \AA^{-1})$")
# plt.ylabel(r"$E(eV)$")
# # plt.title("Add spin term(y direction, j=0.02)")
# # plt.show()
# plt.savefig("no-spin TI film D=6nm.png")
