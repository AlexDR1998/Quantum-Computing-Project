import numpy as np
import math as m
#from gatec import *
from sparse import *

def main():
    # --- Number of qubits and Target Fock value ---
    n = 5
    target = 2
    N = 2**n

    # --- Initialised qubits ---
    q0 = Qubit([1, 0])
    q1 = Qubit([0, 1])

    # --- Initialised gates ---
    I = Identity()
    H = Hadamard(n)   #hadamard all
    X = PauliX(n)   #paulix all
    x = PauliX()
    z = PauliZ()
    cZ = Controlled(z, n)   #controlled z all

    # --- Qreg formation ---
    q = q0&q0
    if n > 2:
        for i in range(n-2):
            q = q0&q

    # --- Number of Iterations calculation ---
    its = int((m.pi/4.0)*(N)**(1/2))


    # --- Fock to Binary Array Conversion ---
    B = [int(x) for x in bin(target)[2:]]
    while len(B) != n:
        B = np.insert(B, 0, 0)
    Binaryform = B

    # --- Oracle PauliX application dependent on Fock ---
    #Initialise
    i = Binaryform[n-1]
    if i == 0:
        Search = x
    elif i == 1:
        Search = I

    #Loop for rest of binary value
    for i in range(n-2, -1, -1):
        if Binaryform[i] == 0:
            Search = x&Search
        elif Binaryform[i] == 1:
            Search = I&Search

    # --- Oracle and Diffusion Gate Calculation ---
    Oracle = Search*cZ*Search
    Diffusion = H*X*cZ*X*H

    # --- Initialising Superposition State ---
    q = H*q

    # --- Grover's Iteration ---
    for i in range(its):
        q = Oracle*q
        q = Diffusion*q

    # --- Measure and Display ---
    q.measure()
    print(q)

main()
