"""
Initial test file for high level implementation of Grover's Algorithm
"""

#NOTE: Values are currently hardcoded; incorporating I/O needed.
#      For now; hard code desired QReg size as n and desired value as desire
#      also hard code number of iterations as its

#NOTE: Coded for understanding, not for efficiency; sectioned so that the steps
#      can be followed easier; for debugging etc.


import numpy as np
import math as m
from gatec import *

def main():

    q0 = Qubit([1, 0])   #|0> qubit
    h = Hadamard()
    I = Identity()
    x = PauliX()

    # --- QReg size ---
    n = 3   #HARDCODE; amend later
    assert type(n) == int, "n must be an integer"
    assert n >= 2, "n must be greater than or equal to 2"

    # --- Formation of QReg ---
    q = q0&q0
    if n > 2:
        for i in range(n-2):
            q = q&q0

    # --- Target ---
    desire = 6   #HARDCODE; amend later
    target = desire - 1  #Fock space value

    # --- Iterations needed ---
    its = 1   #NOTE: This is still a bit buggy when altering number of iterations;
              #      need to go over and check, especially when Grover's iteration
              #      is being sorted.

    # --- Grover gates --- #
    H = Hadamard(n)   #Hadamard all gates
    h1 = I&h   #Hadamard to one gate for n=2
    X = x&x   #PauliX to two gates for n=2
    if n > 2:
        for i in range(n-2):
            h1 = I&h1   #Hadamard to one gate for n>2
            X = x&X   #PauliX to all gates for n>2
    cNOT = CNot()   #for 2 qubit Grover Iteration
    Tof = Toffoli()   #for 3 qubit Grover Iteration
    O = Oracle(n, target)   #Oracle for qubit reg size 2^n and target value in Fock space

    #NOTE: Section to make sure the Qreg gets put to |0>&n

    # --- Formation of superposition from |0>&n state ---
    q = H*q

    # --- Oracle Application ---
    q = O*q

    # --- Grover's Iteration ---
    #NOTE: This section I need to check what parts actually get looped;
    #      I have a funny feeling that its only the bits inside the X's
    if n == 2:
        GI = H * X * h1 * cNOT * h1 * X * H
    elif n == 3:
        GI = H * X * h1 * Tof * h1 * X * H
    q = GI*q

    # --- Measure and Display ---
    q.measure()
    print(q)

main()
