# This program intends to calculate the phase transistion
# of MTI:from C=0 to C=1 when J increase
# Use Dichotomy to Search

from logging import Logger

import numpy as np

import TopoInvCalc as TIC
from TI_Film import Eig, Hamiltonian

# logging.basicConfig(level=logging.INFO)
logger = Logger("DichotomySearch")

settings = {'num_lines': 31,
            'pos_tol':  1e-2,
            'gap_tol': 0.1,
            'move_tol': 0.3,
            'iterator': range(30, 51, 2),
            'min_neighbour_dist': 5e-4,
            }
N = 20
S_ = np.array([[1, 0, 0]]*N)
S = np.array([([s[0]*np.sin(s[1])*np.cos(s[2]), s[0]*np.sin(s[1])
                    * np.sin(s[2]), s[0]*np.cos(s[1])])for s in S_])


def is_Chern(J0, c_tol=0.1):
    h = Hamiltonian(N=N, J=J0, S=S)
    res = TIC.Calc(h, CalcZ2=False)
    return abs(res.Chern) > c_tol
    # Non-trivial Chern Insulator usually C>0.997
    # Trivial Insulator usually C~0.003


def Search(l, r, d_tol=1e-5):
    logger.warning("(l,r)="+str(l)+","+str(r))
    assert l <= r, "Error in the interation: Lower Limit is larger than Upper Limit"

    if (r-l) < d_tol:
        return r
    m = (l+r)/2
    if is_Chern(m):
        return Search(l, m, d_tol=d_tol)
    else:
        return Search(m, r, d_tol=d_tol)


if __name__ == "__main__":
    # print(S_, S)
    res = Search(0, 0.1, d_tol=1e-5)
    with open("Chern PT.txt", "w") as f:
        f.write("Critical J="+str(res))
