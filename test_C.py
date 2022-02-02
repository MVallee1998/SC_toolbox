import ctypes
import numpy.ctypeslib as ctl
from numpy.ctypeslib import ndpointer
import numpy as np

lib = ctypes.cdll.LoadLibrary("./foo.so")
cfoo = lib.cfoo
cfoo.restype = None
cfoo.argtypes = [ctl.ndpointer(np.uint, flags='aligned, c_contiguous')]

pow_2 = np.ones(64,dtype=np.uint)
for k in range(1,64):
    pow_2[k] = pow_2[k-1]*2
cfoo(pow_2)

