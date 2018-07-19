import numpy as np
from scipy import linalg
from matplotlib import pyplot as plt

# Constants related to real material
# unit of length: angstrom
# unit of energy: eV
A1, A2, C, D1, D2 = 2.26, 3.33, -0.0083, 5.74, 30.4
M, B1, B2 = 0.28, 6.86, 44.5
Delta = 3
N = 20  # total layers in z direction, 0.3nm

# Prerequisite matrices
Gzz = np.array(np.diag([D1-B1, D1+B1, D1-B1, D1+B1]), dtype=complex)
Gz = np.zeros([4, 4], dtype=complex)
Gz[0, 1] = Gz[1, 0] = 1.0j*A1
Gz[3, 2] = Gz[2, 3] = -1.0j*A1
t_p = -1*Gzz/(Delta*Delta)-Gz/(2*Delta)
t_m = -1*Gzz/(Delta*Delta)+Gz/(2*Delta)

# print(t_p.T.conjugate()-t_m)
# print(t_m)

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
    hii = E_(kx, ky)*np.identity(4, dtype=complex)+2 * \
        Gzz/(Delta*Delta)+h0  # not include the U(zi)
    # print(hii)

    def row(i):
        return np.column_stack([(hii if j == i else (t_m if j+1 == i else(t_p if j-1 == i else zero)))for j in range(N)])
    return np.row_stack([row(i) for i in range(N)])


# print(Hamiltonian(0.03, 0))
# Simulation parameter, unit = 1/anstrom
# from -xRange to xRange, etc
xRange, yRange, Nx, Ny = 0.1, 0.1, 50, 50
dkx, dky = 2*xRange/Nx, 2*yRange/Ny
bs = np.zeros([Nx+1, Ny+1, 4*N], dtype=float)

# do simulation on ky=0 line
ky = 0
j0 = j = int(Ny/2)
for i in range(Nx+1):
    # kx, ky = -xRange+dkx*i, -yRange+dky*j
    kx = -xRange+dkx*i
    # print(kx, end=", ")
    temp = np.array([x.real for x in (linalg.eig(Hamiltonian(kx, ky))[0])])
    temp.sort()
    # print(temp)
    bs[i, j] = temp
    # print(temp)


X = np.linspace(-xRange, xRange, Nx+1, endpoint=True)
Eig = [([bs[j, j0, i] for j in range(Nx+1)])for i in range(4*N)]


plt.subplot(1, 1, 1)
for y in Eig:
    plt.plot(X, y)
plt.ylim(0, 0.6)
plt.xlim(-0.1, 0.1)
plt.xlabel(r"$k_x(\rm \AA^{-1})$")
plt.ylabel(r"$E(eV)$")

# plt.show()
plt.savefig("Undoped TI film D=6nm.png")

##############################################################
# For debuging
##############################################################

# H = Hamiltonian(0.31, 0)
# print(H)
# la, v = linalg.eig(H)
# print(la)
# Eig = np.array([x.real for x in la])
# print(Eig)
# Eig.sort()
# print(Eig)
# print(Eig[2*N-1])
# print(Eig[2*N])
# H=Hamiltonian(0, 0)

# print(Hamiltonian(0, 0))

# Gz=np.diag([])
# print(Gz)
