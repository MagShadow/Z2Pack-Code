import numpy as np
from scipy import linalg
from matplotlib import pyplot as plt


def Hamiltonian(kx, ky):
    return 0


# Simulation parameter, unit = 1/anstrom
# from -xRange to xRange, etc
N=100
xRange, yRange, Nx, Ny = 0.1, 0.1, 50, 50
dkx, dky = 2*xRange/Nx, 2*yRange/Ny
bs = np.zeros([Nx+1, Ny+1, 4*N], dtype=float)

# do simulation on ky=0 line
ky = 0
j0 = j = int(Ny/2)
for i in range(Nx+1):
    # kx, ky = -xRange+dkx*i, -yRange+dky*j
    kx = -xRange+dkx*i
    print(kx, end="")
    temp = np.array([x.real for x in (linalg.eig(Hamiltonian(kx, ky))[0])])
    temp.sort()
    bs[i, j] = temp
    # print(temp)


X = np.linspace(-xRange, xRange, Nx+1, endpoint=True)
Eig = [([bs[j, j0, i] for j in range(Nx+1)])for i in range(4*N)]


# plt.subplot(1, 1, 1)
# plt.plot(X, Eig[0], color="red")
# plt.plot(X, Eig[1], color="blue")
# plt.plot(X, Eig[2], color="yellow")
# plt.show()

##############################################################
# For debuging
##############################################################
