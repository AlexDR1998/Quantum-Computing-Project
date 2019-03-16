from dense import Qubit as dq
from sparse import Qubit as sq
from lazy import Qubit as lq
from dense import Hadamard as dh
from sparse import Hadamard as sh
from lazy import Hadamard as lh
from dense import Phase as dz
from sparse import Phase as sz
from lazy import Phase as lz
from dense import Controlled as dc
from sparse import Controlled as sc
from lazy import Controlled as lc
from dense import CNot as dx
from sparse import CNot as sx
from lazy import CNot as lx
import numpy as np
import time

from matplotlib import pyplot as plt



def main():
    n = 4
    register_size = []
    for i in range(1,n):
        register_size.append(i)
    densets = []
    for i in range(1,n):
        t1 = time.time()
        #q = dq(i)
        H = dh(i)
        c = dx(i)
        p = dz(np.pi,i)
        q = H&p&c
        
        print(q)
        t2 = time.time()
        timetaken = (t2-t1)
        densets.append(timetaken)
        #print("hi")
    sparsets = []
    for i in range(1,n):
        t1 = time.time()
        #q = sq(i)
        H = sh(i)
        p = sz(np.pi,i)
        c = sx(i)
        q = H&p&c

        print(q)
        t2 = time.time()
        timetaken = (t2-t1)
        sparsets.append(timetaken)
    lazyts = []
    for i in range(1,n):
        t1 = time.time()
        #q = lq(i)
        H = lh(i)
        p = lz(np.pi,i)
        C = lx(i)
        q = H&p&C
        print(q)
        t2 = time.time()
        timetaken = (t2-t1)
        lazyts.append(timetaken)    
   
        
    plt.plot(register_size,densets)
    plt.title("Creating Sparse Registers")
    plt.ylabel("time taken (s)")
    plt.xlabel("size of register")

    
    plt.plot(register_size,sparsets)

    
    plt.plot(register_size,lazyts)
    plt.legend(["dense","sparse","lazy"])
    plt.show()
main()