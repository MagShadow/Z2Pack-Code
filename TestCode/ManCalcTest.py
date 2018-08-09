import numpy as np
import os

C = []
with open("wcc_result.txt", "r") as f:
    for line in f.readlines():
        C = C+str.split(line.strip())

C = np.array([float(x) for x in C[:int(len(C)/2)]])
print(C)
print(len(C))
print(sum(C))
