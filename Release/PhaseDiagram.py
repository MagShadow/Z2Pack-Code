# This define a common phase diagram with 2 axes

import numpy as np
# 在服务器端没有DISPLAY资源，因此必须换成无交互的后端。默认为"TkAgg"。
import matplotlib as mpl
mpl.use("Agg")
from matplotlib import pyplot as plt

from copy import deepcopy
from json import dumps, loads
from itertools import product
from datetime import datetime
from multiprocessing import Manager

from _utils import nt, settings, settings_strict


class PhaseDiag(object):
    r'''
    Phase Diagram for different parameters.

    Func `run` accept a function, two axes to draw the phase diagram.
    '''

    def __init__(self, PD=None):
        if PD == None:
            return
        assert isinstance(
            PD, PhaseDiag), "Can only be initialized from a PhaseDiag object!"
        self.func = PD.func
        self.data = PD.data
        self.info = PD.info

    def run(self, func, X_L, X_H, X_N, Y_L, Y_H, Y_N, X_Name="", Y_Name=""):
        self.T_start = datetime.now()
        print("Start Calculation at ", str(self.T_start))
        self.func = func
        X = np.linspace(X_L, X_H, num=X_N, endpoint=True)
        Y = np.linspace(Y_L, Y_H, num=Y_N, endpoint=True)

        m = Manager()
        p = m.Pool()
        Phase = m.list([m.list([0]*Y_N) for i in range(X_N)])

        for i, j in product(list(range(X_N)), list(range(Y_N))):
            _X, _Y = X[i], Y[j]
            p.apply_async(self.func, args=(_X, _Y, i, j, Phase))

        p.close()
        p.join()

        self.data = [list(x) for x in list(Phase)]
        self.info = dict(X_L=X_L, X_H=X_H, X_N=X_N, Y_L=Y_L,
                         Y_H=Y_H, Y_N=Y_N, X_Name=X_Name, Y_Name=Y_Name)

        self.T_end = datetime.now()
        print("End Calculation at ", str(self.T_end))
        print("Total time:", str(self.T_end-self.T_start))

        return self

    def read(self, filename=""):
        assert filename != "", "filename cannot be empty!"
        with open(filename, "r") as f:
            d = f.read()
            self.info = loads(d)
            self.data = deepcopy(self.info["data"])
            del self.info["data"]

    def write(self, filename="", mode="w"):
        assert mode in ["w", "a"], "Writing Mode can only be 'w' or 'a'!"
        if filename == "":
            filename = "PhaseDiag_Data_"+nt()+".txt"
        elif ((len(filename) <= 4) or (filename[-4:] != ".txt")):
            filename = filename+".txt"
        d = deepcopy(self.info)
        d["data"] = self.data
        with open(filename, mode) as f:
            f.write(dumps(d))

    def draw(self, title="", xlabel="", ylabel="", filename=""):
        X_L, X_H, X_N = self.info["X_L"], self.info["X_H"], self.info["X_N"]
        Y_L, Y_H, Y_N = self.info["Y_L"], self.info["Y_H"], self.info["Y_N"]

        X = np.linspace(X_L, X_H, num=X_N+1, endpoint=True)
        Y = np.linspace(Y_L, Y_H, num=Y_N+1, endpoint=True)
        Z = self.data

        fig, ax = plt.subplots()
        cmap = mpl.colors.ListedColormap(["r", "g", "b", "c"])
        norm = mpl.colors.BoundaryNorm(list(range(5)), cmap.N)
        x, y = np.meshgrid(X, Y, indexing="ij")
        ax.pcolormesh(x, y, np.array(Z), cmap=cmap, vmin=0, vmax=3)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax2 = fig.add_axes([0.92, 0.1, 0.03, 0.8])
        mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm)

        if filename != "":
            plt.savefig(filename+".png")
        else:
            plt.show()
