import numpy as np
import os
import time
import random
from datetime import datetime
from multiprocessing import Pool

from TI_Film import Hamiltonian
import TopoInvCalc as TIC

random.seed()


def RdSch(index):
    N, J = 20, 0.02
    S_ = np.zeros([N, 3])
    for i in range(N):
        S_[i, 0] = 1
        S_[i, 1] = random.random()*np.pi
        S_[i, 2] = random.random()*np.pi*2
    if index==0:
        S_=np.array([[1,0,0]]*N)
    S = np.array([([s[0]*np.sin(s[1])*np.cos(s[2]), s[0]*np.sin(s[1])
                    * np.sin(s[2]), s[0]*np.cos(s[1])])for s in S_])
    h = Hamiltonian(N=N, J=J, S=S)
    res=TIC.Calc(h, CalcZ2=True)
    if (abs(res.Chern)>0.1) or res.Z2:
        p=os.path.join("SearchResult","Result-"+str(index)+" "+str(datetime.now()))
        with open(p,"w") as f:
            f.write("Chern="+res.Chern+"\n")
            f.write("Z2="+res._Z2+"\n")
            f.write("Spin Distribution:",str(S))



# def Multi_Test(name):
#     print('Run task %s (%s)...' % (name, os.getpid()))
#     start = time.time()
#     time.sleep(random.random() * 3)
#     end = time.time()
#     print('Task %s runs %0.2f seconds.' % (name, (end - start)))


def Search(N):
    p = Pool()
    for i in range(N):
        p.apply_async(RdSch, args=(i,))
    p.close()
    p.join()
    print("Done!")

def mkdir(path):
    path=path.strip()
    path=path.rstrip("\\")
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)
    return

if __name__ == "__main__":
    mkdir("SearchResult")
    Search(100)
    # print('Parent process %s.' % os.getpid())
    # p = Pool()
    # for i in range(20):
    #     p.apply_async(Multi_Test, args=(i,))
    # print('Waiting for all subprocesses done...')
    # p.close()
    # p.join()
    # print('All subprocesses done.')
