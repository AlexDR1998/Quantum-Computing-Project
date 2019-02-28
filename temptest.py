from lazy import *
#from sparse import *
#from gatec import *
import numpy as np
import time


def main():
    h = Hadamard()
    t1 = time.time()
    H = h&h&h&h&h&h&h&h&h
    t2 = time.time()
    print(H)
    print(t2-t1)
main()