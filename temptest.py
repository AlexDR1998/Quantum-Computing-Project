from gatec import QMatrix as d
from sparse import SMatrix as s
from lazy import LMatrix as l
import numpy as np
import time
from matplotlib import pyplot as plt


def main():
    n = 7
    register_size = []
    for i in range(n):
        register_size.append(i)
    densets = []
    for i in range(n):
        t1 = time.time()
        h = d.Hadamard(i)
        print(h)
        t2 = time.time()
        time = (t2-t1)
        densets.append(time)
    sparsets = []
    for i in range(n):
        t1 = time.time()
        h = s.Hadamard(i)
        print(h)
        t2 = time.time()
        time = (t2-t1)
        sparsets.append(time)
    lazyts = []
    for i in range(n):
        t1 = time.time()
        h = l.Hadamard(i)
        print(h)
        t2 = time.time()
        time = (t2-t1)
        lazyts.append(time)
        
    plt.plot(register_size,densets)
    plt.title("dense implementation")
    plt.ylabel("time to make hadamard")
    plt.xlabel("size of qubit register acted on")
    plt.show()
    
    plt.plot(register_size,sparsets)
    plt.title("dense implementation")
    plt.ylabel("time to make hadamard")
    plt.xlabel("size of qubit register acted on")
    plt.show()
    
    plt.plot(register_size,lazyts)
    plt.title("dense implementation")
    plt.ylabel("time to make hadamard")
    plt.xlabel("size of qubit register acted on")
    plt.show()
main()