import numpy as np
import os
import time
import random
from json import dumps, loads
from sys import stderr
from datetime import datetime
from multiprocessing import Pool, Queue, Manager

from TI_Film import Hamiltonian
import TopoInvCalc as TIC

random.seed()


def nt():
    return datetime.now().strftime("%y-%m-%d-%H-%M-%S")


sfname = os.path.join("SearchResult", nt())


def RdSch(index, q=None):
    random.seed()

    print("Index %d started at %s" % (index, datetime.now()))
    N, J = 20, 0.02
    S_ = np.zeros([N, 3])
    for i in range(N):
        S_[i, 0] = 1
        S_[i, 1] = random.random()*np.pi
        S_[i, 2] = random.random()*np.pi*2
    if index == 0:
        S_ = np.array([[1, 0, 0]]*N)
    S = np.array([([s[0]*np.sin(s[1])*np.cos(s[2]), s[0]*np.sin(s[1])
                    * np.sin(s[2]), s[0]*np.cos(s[1])])for s in S_])
    S = S*random.random()
    # print(S)
    h = Hamiltonian(N=N, J=J)
    res = TIC.Calc(h, CalcZ2=True, LogOut=False)
    # print(res.Chern)
    if res.Z2:
        # print("True!")
        d = dict(index=index, C=res.Chern, Z2=res._Z2, SpinDist=S.tolist())
        if q != None:
            # print("Push!")
            q.put(dumps(d))
        else:
            mkdir(sfname)
            with open(os.path.join(sfname, "Result-" + str(index)+".txt"), "w") as f:
                f.write(d)
                # q.push("Index="+str(i)+" ,Chern=" +
                #        str(res.Chern)+" ,Z2="+str(res._Z2)+"\n")
                # p = os.path.join("SearchResult", "Result-" +
                #                  str(index)+".txt")
                # with open(p, "w") as f:
                #     f.write("Chern="+str(res.Chern)+"\n")
                #     f.write("Z2="+str(res._Z2)+"\n")
                #     f.write("Spin Distribution:"+str(S))

                # def Multi_Test(name):
                #     print('Run task %s (%s)...' % (name, os.getpid()))
                #     start = time.time() #     time.sleep(random.random() * 3)
                #     end = time.time()
    print("Index %d ended at %s, Chern=%.3f, Z2=%s" %
          (index, datetime.now(), res.Chern, str(res._Z2)))


def Search(N):
    p, m = Pool(), Manager()
    q = m.Queue()
    for i in range(N):
        p.apply_async(RdSch, args=(i, q,))
    p.close()
    p.join()
    res = []
    # print(q.empty())
    while not q.empty():
        res.append(q.get())
    if res != []:
        r_index = []
        with open(sfname+".txt", "w") as f:
            for r in res:
                f.write(r+"\n")
                # print(loads(r))
                r_index.append(loads(r).get("index"))
        with open(os.path.join("SearchResult", "Result.txt"), "a") as f:
            f.write(sfname+" "+str(r_index)+"\n")
        # print(r_index)
    else:
        print("No Result!")
        with open(os.path.join("SearchResult", "Result.txt"), "a") as f:
            f.write(sfname+" "+str([])+"\n")
    stderr.write("Done!"+nt())


def mkdir(path):
    path = path.strip()
    path = path.rstrip("\\")
    isExists = os.path.exists(path)
    if not isExists:
        os.makedirs(path)
    return


if __name__ == "__main__":
    Search(300)
    # print('Parent process %s.' % os.getpid())
    # p = Pool()
    # for i in range(20):
    #     p.apply_async(Multi_Test, args=(i,))
    # print('Waiting for all subprocesses done...')
    # p.close()
    # p.join()
    # print('All subprocesses done.')
