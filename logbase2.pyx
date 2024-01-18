from math import log2
cpdef float logbase2(float x):
    if x == 0:
        return 0
    else:
        return log2(x)