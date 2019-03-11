from dense import Hadamard as d
from sparse import Hadamard as s
from lazy import Hadamard as l
import numpy as np
import time
from matplotlib import pyplot as plt


def main():
    n = 10
    
    register_size = []
    for i in range(n):
        register_size.append(i)
    densets = []
    for i in range(n):
        t1 = time.time()
        h = d(i)
        print(h)
        t2 = time.time()
        timetaken = (t2-t1)
        densets.append(timetaken)
    sparsets = []
    for i in range(n):
        t1 = time.time()
        h = s(i)
        print(h)
        t2 = time.time()
        timetaken = (t2-t1)
        sparsets.append(timetaken)
    lazyts = []
    for i in range(n):
        t1 = time.time()
        h = l(i)
        print(h)
        t2 = time.time()
        timetaken = (t2-t1)
        lazyts.append(timetaken)
        
    plt.plot(register_size,densets)
    plt.title("dense implementation")
    plt.ylabel("time to make hadamard")
    plt.xlabel("size of qubit register acted on")
    plt.show()
    
    plt.plot(register_size,sparsets)
    plt.title("sparse implementation")
    plt.ylabel("time to make hadamard")
    plt.xlabel("size of qubit register acted on")
    plt.show()
    
    plt.plot(register_size,lazyts)
    plt.title("lazy implementation")
    plt.ylabel("time to make hadamard")
    plt.xlabel("size of qubit register acted on")
    plt.show()
main()