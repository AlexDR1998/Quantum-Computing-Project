import numpy as np
import math as m
from gatec import *

def main():

    q1 = Qubit([1, 0])
    q2 = Qubit([0, 1])
    I = Identity()

    n = 4
    N = 2**n
    qreg = q1&q1&q1&q1
    H = Hadamard(4)
    X = PauliX(4)
    x = PauliX(1)

    its = int((m.pi/4.0)*(N)**(1/2))

    Binaryform = [int(x) for x in bin(n)[2:]]

    i = Binaryform[n-1]
    if i = 0:
        Search = x
    elif i = 1:
        Search = I

    for i in range(n-2, 0, -1):
        if i = 0:
            Search = x&O
        elif i = 1:
            Search = I&O

    qf = H*qreg

    for i in range(its):
        


    qf.measure()
    print(qf)

main()
