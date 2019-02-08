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

    # --- HARDCODED VALUES; Adjust to I/O later ---
    n = 3
    target = 4

    # --- QReg size ---
    N = 2**n
    assert type(n) == int, "n must be an integer"
    assert n >= 2, "n must be greater than or equal to 2"

    # --- Formation of QReg ---
    q = q0&q0
    if n > 2:
        for i in range(n-2):
            q = q0&q

    # --- Target ---
    TAR = target - 1  #Fock space value is 1 above indice value

    # --- Iterations needed ---
    its = int((m.pi/4.0)*(N)**(1/2))
    #NOTE: This is still a bit buggy; the target isn't achieved everytime.
    #      Need to go over and check, especially when Grover's iteration is being sorted.

    # --- Grover gates --- #
    H = Hadamard(n)   #Hadamard all gates
    h1 = I&h   #Hadamard to one gate for n=2
    if n > 2:
        for i in range(n-2):
            h1 = I&h1   #Hadamard to one gate for n>2
    X = PauliX(n)   #PauliX to n qubits
    cNOT = CNot()   #for 2 qubit Grover Iteration
    Tof = Toffoli()   #for 3 qubit Grover Iteration
    O = Oracle(n, TAR)   #Oracle for qubit reg size 2^n and target value in Fock space

    #NOTE: Section to make sure the Qreg gets put to |0>^&n

    # --- Formation of superposition from |0>&n state ---
    q = H*q

    # --- Grover's Iteration ---
    #NOTE: This section I need to check what parts actually get looped

    for i in range(its):
        q = O*q   #Oracle application
        if n == 2:
            GI = H * X * h1 * cNOT * h1 * X * H
        elif n == 3:
            GI = H * X * h1 * Tof * h1 * X * H
        q = GI*q

    # --- Measure and Display ---
    q.measure()
    print(q)

main()
