# Define some useful constants and functions

from datetime import datetime
import numpy as np
settings = {'num_lines': 51,
            'pos_tol':  1e-2,
            'gap_tol': 0.05,
            'move_tol': 0.2,
            'iterator': range(50, 81, 4),
            'min_neighbour_dist': 1e-4,
            }
settings_strict = {'num_lines': 81,
                   'pos_tol':  1e-3,
                   'gap_tol': 0.001,
                   'move_tol': 0.4,
                   'iterator': range(80, 201, 5),
                   'min_neighbour_dist': 5e-6,
                   }

# sig_xyz = pauli_xyz tensor identity
sig_x, sig_y, sig_z = [np.zeros([4, 4], dtype=complex) for i in range(3)]
sig_x[0, 2], sig_x[1, 3], sig_x[2, 0], sig_x[3, 1] = 1, 1, 1, 1
sig_y[0, 2], sig_y[1, 3], sig_y[2, 0], sig_y[3, 1] = -1j, -1j, 1j, 1j
sig_z[0, 0], sig_z[1, 1], sig_z[2, 2], sig_z[3, 3] = 1, 1, -1, -1
sig = np.array([sig_x, sig_y, sig_z])

# Constants related to real material
# unit of length: angstrom
# unit of energy: eV

# Consts from PHYSICAL REVIEW B 82, 045122 (2010)
CONST_HJZ = {
    "A1": 2.26, "A2": 3.33, "C": -0.0083, "D1": 5.74, "D2": 30.4,
    "M": 0.28, "B1": 6.86, "B2": 44.5
}
CONST_HZL_O = {  # Interesting!!!!!
    "A1": 3.3, "A2": 4.1, "C": -0.0068, "D1": 1.2, "D2": -30.1,
    "M": 0.28, "B1": 1.5, "B2": -54.1
}
# Consts from New J. Phys. 12 043048 (2010)
# B1 and D1 has different signs with their papers, because the M term and E term
# in different papers have different forms
CONST_HZL = {
    "A1": 3.3, "A2": 4.1, "C": -0.0068, "D1": -1.2, "D2": -30.1,
    "M": 0.28, "B1": -1.5, "B2": -54.1
}
# Consts from New J. Phys. 12 043048 (2010)
# B1 and D1 has different signs with their papers, because the M term and E term
# in different papers have different forms
CONST_HZL_M = {  # Modified for Some Testing.
    "A1": 3.3, "A2": 4.1, "C": -0.0068, "D1": 1.2, "D2": 30.1,
    "M": 0.28, "B1": 1.5, "B2": 54.1
}


def nt():
    '''
    Format the time.
    '''
    return datetime.now().strftime("%y-%m-%d-%H-%M-%S")


def convertSpin(S):
    r'''
    Convert the spin array with a form of $(r,\theta,\phi)$ into a form of $(Sx,Sy,Sz)$
    '''
    return np.array([([s[0]*np.sin(s[1])*np.cos(s[2]), s[0]*np.sin(s[1])
                           * np.sin(s[2]), s[0]*np.cos(s[1])])for s in S])
