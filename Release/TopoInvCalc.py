import z2pack
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
from matplotlib import pyplot as plt
import os
import logging
from datetime import datetime

from TI_Film import Eig, plotBS, Hamiltonian


# Constant for Z2Pack Calculation
identity = np.identity(2, dtype=complex)
pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
pauli_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
pauli_z = np.array([[1, 0], [0, -1]], dtype=complex)
pauli_vector = list([pauli_x, pauli_y, pauli_z])

settings = {'num_lines': 51,
            'pos_tol':  1e-4,
            'gap_tol': 0.001,
            'move_tol': 0.5,
            'iterator': range(50, 101, 5),
            'min_neighbour_dist': 5e-6,
            }


def TopoOrder(res, _Chern_tol=0.1):
    '''
    Accept an Result from TopoInvCalc,
    return
        0 if trivial;
        1 if Z2 but not Chern;
        2 if Chern but not Z2;
        3 if Chern and Z2;
    '''
    C, Z = 1 if abs(res.Chern) > _Chern_tol else 0, int(res.Z2)
    return C*2+Z


def Calc(ham, bands=None, surf=lambda k1, k2: [k1-0.5, k2-0.5], KScale=1, CalcZ2=True, LogOut=True, Timer=True, settings=settings):
    if not LogOut:
        logging.getLogger('z2pack').setLevel(logging.ERROR)

    def h0(k): return ham(k[0]/KScale, k[1]/KScale)

    class Res(object):

        def __init__(self):
            if Timer:
                T_start = datetime.now()

            s0 = z2pack.hm.System(h0, bands=bands, dim=2)
            result = z2pack.surface.run(
                system=s0,
                # parameter of surface is moduled by 2pi
                surface=surf,
                # surface=lambda k1, k2: [k2, k1/2],
                **settings
                # save_file="savefile.msgpack"
            )
            self.result = result
            self.Chern = z2pack.invariant.chern(result)
            self.CalcZ2 = CalcZ2
            if CalcZ2 == True:
                self.result_l = z2pack.surface.run(
                    system=s0, surface=lambda k1, k2: surf(k1/2, k2), **settings)
                self.result_r = z2pack.surface.run(
                    system=s0, surface=lambda k1, k2: surf(k1/2+0.5, k2), **settings)
                Z2_1, Z2_2, Z2_T = z2pack.invariant.z2(self.result_l), z2pack.invariant.z2(
                    self.result_r), z2pack.invariant.z2(self.result)
                self.Z2 = True if (Z2_1*Z2_2*(1-Z2_T) == 1) else False
                self._Z2 = [Z2_1, Z2_2, Z2_T]

            if Timer:
                T_end = datetime.now()
                print("Running Time:"+str(T_end-T_start))

        def plotChern(self, title="", filename="", start=0, end=1):
            fig, ax = plt.subplots()
            z2pack.plot.wcc(self.result, axis=ax)
            ax.set_xlabel(r"$k_x$")
            ax.set_ylabel(r"$\bar{x}$")
            ax.set_xlim(start, end)
            ax.set_title("Chern="+str(self.Chern))
            fig.suptitle(title)
            if filename != "":
                path = os.path.join("Pictures", filename+".png")
                plt.savefig(path)
            else:
                plt.show()

        def plotZ2(self, title="", filename="", start=0, end=1, start_l=0, end_l=1, start_r=0, end_r=1):
            if self.CalcZ2 != True:
                raise Exception("Haven't Calculated Z2")
            else:
                fig, ax = plt.subplots(1, 3)
                z2pack.plot.wcc(self.result_l, axis=ax[0])
                z2pack.plot.wcc(self.result_r, axis=ax[1])
                z2pack.plot.wcc(self.result, axis=ax[2])
                for i in range(3):
                    ax[i].set_xlabel(r"$k_x$")
                    ax[i].set_ylabel(r"$\bar{x}$")
                    ax[i].set_title("Z2="+str(self._Z2[i]))

                ax[0].set_xlim(start_l, end_l)
                ax[1].set_xlim(start_r, end_r)
                ax[2].set_xlim(start, end)

                fig.suptitle(title)
                if filename != "":
                    path = os.path.join("Pictures", filename+".png")
                    plt.savefig(path)
                else:
                    plt.show()

    return(Res())


def Calc_Man(ham, bands=None, surf=lambda k1, k2: [k1-0.5, k2-0.5], KScale=1, CalcZ2=True):
    '''
    Try to calculate the Z2 index manually by calculating the Chern Number of Each Band.
    Still under constructing
    '''
    def h0(k): return ham(k[0]/KScale, k[1]/KScale)
    s0 = z2pack.hm.System(h0, bands=bands, dim=2)
    result = z2pack.surface.run(
        system=s0,
        # parameter of surface is moduled by 2pi
        surface=surf,
        # surface=lambda k1, k2: [k2, k1/2],
        **settings
        # save_file="savefile.msgpack"
    )
    X, _wcc = result.t, result.wcc
    L, N = len(_wcc), len(_wcc[0])
    wcc = np.array([[_wcc[j][i] for j in range(L)]for i in range(N)])

    logger = logging.Logger("BandCalc")
    res = [0]*N
    for i in range(N):
        res[i] = z2pack.surface.run(system=z2pack.hm.System(
            h0, dim=2, bands=[i]), surface=surf, **settings)
        logger.warning("Band="+str(i))

    c = np.array([z2pack.invariant.chern(r) for r in res])
    print(c)

    print("Num of WCC~0.5 at k=0.5", sum(
        [(1 if abs(w[int((L-1)/2)]-0.5) < 1e-7 else 0)for w in wcc]))
    print("C, sum of C", z2pack.invariant.chern(result), sum(c))
    fig, ax = plt.subplots()
    z2pack.plot.wcc(result, axis=ax)
    plt.savefig("TestWCC_Alone_Z-Spin.png")
