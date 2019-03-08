from lazy import *
#from lazzz import *
#from sparse import *
#from gatec import *
import numpy as np
import time


def main():
    h = Hadamard(2)
    q = Qubit(2)
    X = PauliX()
    
    print(h)
    print(q)
    t1 = time.time()
    #gate = 
    a = h*q
    
    print(a)
    t2 = time.time()
    print(t2-t1)
main()