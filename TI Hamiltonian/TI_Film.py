# Define Hamiltonian for TI Film

import numpy as np
from scipy import linalg
from matplotlib import pyplot as plt
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
N_ = 20  # total layers in z direction, N=20->d=0.6nm
J_ = 0.02  # spin coupling energy

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


# vS[i] discribe the spin distribution
# first component is the intensity, from 0 to 1
# second and three term is the theta and phi, describe the direction of the spin
# S_[i] is the spin vector
vS = np.zeros([N_, 3])
for i in range(N_):
    vS[i, 0] = 1
    vS[i, 1] = 0
    vS[i, 2] = 0

S_ = np.array([([s[0]*np.sin(s[1])*np.cos(s[2]), s[0]*np.sin(s[1])
                 * np.sin(s[2]), s[0]*np.cos(s[1])])for s in vS])
U_ = np.zeros([N_])

xRange, yRange, Nx, Ny = 0.05, 0.05, 50, 50
X_ = np.linspace(-xRange, xRange, Nx+1, endpoint=True)
Y_ = np.linspace(-yRange, yRange, Ny+1, endpoint=True)

# Return a 4N*4N matrix
def Hamiltonian(N=N_, J=J_, S=S_, U=U_):
    def _ham(kx, ky):
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

        h_ = np.array([((E_(kx, ky)+U[i])*np.identity(4, dtype=complex)+2 *
                        Gzz/(Delta*Delta)+h0-J*sum([S[i, j]*sig[j] for j in range(3)])) for i in range(N)])

        def row(i):
            return np.column_stack([(h_[i] if j == i else (t_m if j+1 == i else(t_p if j-1 == i else zero)))for j in range(N)])
        return np.row_stack([row(i) for i in range(N)])
    return _ham


def Eig(h, xRange=xRange, yRange=yRange, Nx=Nx, Ny=Ny):
    '''
        Accept callable object h(kx,ky)
    '''
    N1, N2 = h(0, 0).shape
    N = int(N1/4)
    dkx, dky = 2*xRange/Nx, 2*yRange/Ny
    bs = np.zeros([Nx+1, Ny+1, 4*N], dtype=float)
    # bs = [[0]*(Ny+1)]*(Nx+1)
    for i in range(Nx+1):
        for j in range(Ny+1):
            kx, ky = -xRange+dkx*i, -yRange+dky*j
            temp = np.array([x.real for x in (linalg.eig(h(kx, ky))[0])])
            temp.sort()
            bs[i, j] = temp

    return(np.array([([([bs[i, j, n] for j in range(Ny+1)]) for i in range(Nx+1)])
                     for n in range(4*N)]))





def plotBS(E, start, end, X=X_, Y=Y_, filename="",title=""):
    x, y = np.meshgrid(X, Y)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for Z in E[start:end]:
        # print(Z)
        ax.plot_surface(x, y, Z, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)

    # fig.colorbar(surf, shrink=0.5, aspect=5)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    ax.set_xlabel(r"$k_x(\rm \AA^{-1})$")
    ax.set_ylabel(r"$k_y(\rm \AA^{-1})$")
    ax.set_zlabel(r"$E/\rm{eV}$")
    ax.set_title(title)
    if filename == "":
        plt.show()
    else:
        path = os.path.join("Pictures", filename+".png")
        plt.savefig(path)


if __name__ == "__main__":
    N, J = 20, 0.02
    S_ = np.zeros([N, 3])
    for i in range(N):
        S_[i, 0] = 1
        S_[i, 1] = 0
        S_[i, 2] = 0
    S = np.array([([s[0]*np.sin(s[1])*np.cos(s[2]), s[0]*np.sin(s[1])
                    * np.sin(s[2]), s[0]*np.cos(s[1])])for s in S_])
    h = Hamiltonian(N=N, J=J, S=S)
    e = Eig(h)

    plotBS(e, 2*N-2, 2*N+2,title="TI Film: No Spin Term, J=0.02, 4 bands")
