# Define Hamiltonian for TI Film
import os
import numpy as np
from scipy import linalg
from numbers import Integral
import matplotlib as mpl
mpl.use("Agg")
from matplotlib import pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from _utils import sig, convertSpin


def SpinZ(M):
    r'''
    Return an array of Spin in Z direction, with a form of $(Sx,Sy,Sz)$
    '''
    # vS[i] discribe the spin distribution
    # first component is the intensity, from 0 to 1
    # second and three term is the theta and phi, describe the direction of the spin
    # S_[i] is the spin vector
    vS = np.array([[1, 0, 0]]*M)
    S_ = np.array([([s[0]*np.sin(s[1])*np.cos(s[2]), s[0]*np.sin(s[1])
                     * np.sin(s[2]), s[0]*np.cos(s[1])])for s in vS])
    return S_


def Hamiltonian(N=20, J=0, S=[], U=[], Delta=3,
                *, A1=2.26, A2=3.33, C=-0.0083, D1=5.74, D2=30.4,
                M=0.28, B1=6.86, B2=44.5):
    '''
    Default constants from PHYSICAL REVIEW B 82, 045122 (2010)
    '''
    assert isinstance(N, Integral), "N should be an interger!"
    assert (J > 0) or (abs(J) < 1e-9), "J should >0!"
    if S == []:
        S = SpinZ(N)
    if U == []:
        U = np.zeros([N])
    assert len(S) == N, "Length of S distribution should equal to N"
    assert len(U) == N, "Length of U distribution should equal to N"

    # Prerequisite matrices
    Gzz = np.array(np.diag([D1-B1, D1+B1, D1-B1, D1+B1]), dtype=complex)
    Gz = np.zeros([4, 4], dtype=complex)
    Gz[0, 1] = Gz[1, 0] = 1.0j*A1
    Gz[3, 2] = Gz[2, 3] = -1.0j*A1
    t_p = -1*Gzz/(Delta*Delta)-Gz/(2*Delta)
    t_m = -1*Gzz/(Delta*Delta)+Gz/(2*Delta)

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


xRange, yRange, Nx, Ny = 0.05, 0.05, 50, 50
X_ = np.linspace(-xRange, xRange, Nx+1, endpoint=True)
Y_ = np.linspace(-yRange, yRange, Ny+1, endpoint=True)


def Eig(h, xRange=xRange, yRange=yRange, Nx=Nx, Ny=Ny):
    '''
        Accept callable object h(kx,ky);
        Return a matrix, Eig[Band_Index,kx_index,ky_index]
    '''
    assert isinstance(Nx, Integral), "Nx should be an interger!"
    assert isinstance(Ny, Integral), "Ny should be an interger!"

    N1 = h(0, 0).shape[0]
    N = int(N1/4)
    dkx, dky = 2*xRange/Nx, 2*yRange/Ny
    bs = np.zeros([Nx+1, Ny+1, 4*N], dtype=float)
    for i in range(Nx+1):
        for j in range(Ny+1):
            kx, ky = -xRange+dkx*i, -yRange+dky*j
            temp = np.array([x.real for x in (linalg.eig(h(kx, ky))[0])])
            temp.sort()
            bs[i, j] = temp

    temp = [([([bs[i, j, n] for j in range(Ny+1)]) for i in range(Nx+1)])
            for n in range(4*N)]
    return np.array(temp)


def Gap(E, k=None):
    if k != None:
        assert type(k) == tuple, "Accept para like (kx_index,ky_index)!"
        assert len(k) == 2, "Accept para like (kx_index,ky_index)!"
        N = int(E.shape[0]/4)
        _Gap = E[2*N]-E[2*N-1]
        return min(_Gap)
    else:
        return E[k[0], k[1]]


def plotBS(E, start, end, X=X_, Y=Y_, filename="", title=""):
    x, y = np.meshgrid(X, Y, indexing="ij")
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for Z in E[start:end]:
        ax.plot_surface(x, y, Z, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)

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


def plotLine(h, start, end, xRange=0.05, Nx=20, axis="x", filename="", title=""):
    assert isinstance(Nx, Integral), "N should be an integer!"

    step = 2*xRange/Nx
    if axis == "x":
        X = np.array([[-xRange+step*i, 0]for i in range(Nx+1)])
    else:
        X = np.array([[0, -xRange+step*i]for i in range(Nx+1)])
    N_band = h(0, 0).shape[0]
    E = np.zeros([Nx+1, N_band], dtype=float)
    for i in range(Nx+1):
        temp = np.array([x.real for x in (linalg.eig(h(X[i, 0], X[i, 1]))[0])])
        temp.sort()
        E[i] = temp
    # print(Eig)
    Z = np.array([[E[i, j] for i in range(Nx+1)] for j in range(N_band)])
    # print(Z)
    plt.subplot(1, 1, 1)
    x = np.linspace(-xRange, xRange, Nx+1, endpoint=True)

    for b in Z[start:end]:
        plt.plot(x, b)
    plt.xlabel(r"$k_x(\rm \AA^{-1})$")
    plt.ylabel(r"$E(eV)$")
    plt.title(title)
    if filename != "":
        plt.savefig(filename+".png")
    else:
        plt.show()
    # print("End!")


if __name__ == "__main__":
    N, Delta, J = 18, 3.33, 0.00

    # _S = np.zeros([N, 3])
    # S = convertSpin(_S)

    xRange, yRange, Nx, Ny = 0.05, 0.05, 30, 30
    X_ = np.linspace(-xRange, xRange, Nx+1, endpoint=True)
    Y_ = np.linspace(-yRange, yRange, Ny+1, endpoint=True)

    h = Hamiltonian(N=N, J=J, Delta=Delta)
    plotLine(h, 2*N-2, 2*N+2, xRange=xRange, Nx=Nx)
    # e = Eig(h, xRange=xRange, yRange=yRange, Nx=Nx, Ny=Ny)
    # plotBS(e, 2*N-2, 2*N+2, X=X_, Y=Y_,
    #        title="TI Film, Spin: +z, J=0.02, 4 bands")
