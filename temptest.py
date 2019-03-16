from dense import Qubit as dq
from sparse import Qubit as sq
from lazy import Qubit as lq
from dense import Hadamard as dh
from sparse import Hadamard as sh
from lazy import Hadamard as lh
from dense import Phase as dz
from sparse import Phase as sz
from lazy import Phase as lz
from dense import CNot as dc
from sparse import CNot as sc
from lazy import CNot as lc
from dense import PauliX as dx
from sparse import PauliX as sx
from lazy import PauliX as lx
import numpy as np
import time
import sys


from matplotlib import pyplot as plt



def main():
    n = 8
    register_size = []
    for i in range(1,n):
        register_size.append(i)
    densets = []
    for i in range(1,n):
        t1 = time.time()
        h = dc(i)
        

      
        
        print(h)
        t2 = time.time()
        timetaken = t2-t1
        densets.append(timetaken)
        #print("hi")
    sparsets = []
    for i in range(1,n):
        t1 = time.time()
        h = sc(i)
        
        print(h)
        t2 = time.time()
        time_taken = t2-t1
        sparsets.append(timetaken)
    lazyts = []
    for i in range(1,n):
        t1 = time.time()
        h = lc(i)

     
        print(h)
        t2 = time.time()
        timetaken = t2-t1
        lazyts.append(timetaken)    
   
        
    plt.plot(register_size,densets)
    plt.title("Tensoring large gates together")
    plt.ylabel("time taken (s)")
    plt.xlabel("size of gates")

    print(sparsets)
    plt.plot(register_size,sparsets)

    
    plt.plot(register_size,lazyts)
    plt.legend(["dense","sparse","lazy"])
    plt.show()
main()