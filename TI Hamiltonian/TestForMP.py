import random
import time
from functools import partial
from multiprocessing import Lock, Pool, Manager


def Run(i, Y, lock):
    print("Run: index=", i)
    time.sleep(int(random.random()*10))
    # lock.acquire()
    # print("Get Lock!, index=", i)
    Y[i] = i*i
    print("Write: index=", i)
    # lock.release()


def Run_without_lock(i, Y):
    print("Run: index=", i)
    # time.sleep(int(random.random()*10))
    # lock.acquire()
    # print("Get Lock!, index=", i)
    Y[i] = i*i
    print("Write: index=", i)


if __name__ == "__main__":
    p, m = Pool(), Manager()
    # l = m.Lock()
    # global Y
    Y = m.list([0]*10)
    # pRun = partial(Run, lock=l)
    for i in range(10):
        p.apply_async(Run_without_lock, args=(i, Y,))
    p.close()
    p.join()
    print(Y)
