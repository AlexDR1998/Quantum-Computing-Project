import numpy as np
import math as m
import time as t
from gatec import *
#from sparse import *

def main():
    # --- Number of qubits and Target Fock value ---
    n = int(input('How many qubits? '))
    assert type(n) == int, "n must be an integer greater than or equal to 2"
    assert n >= 2, "n must be an integer greater than or equal to 2"
    N = 2**n
    target = int(input('What Fock space value would you like to find? '))
    assert type(target) == int, "Target must be an integer greater than or equal to 2"
    assert target >= 0, "Target must be an integer greater than or equal to 0"
    assert target <= N-1, "Target must be an integer less than or equal to " + str(N-1)

    # --- Timer Initialise ---
    init = t.time()

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

    # --- Oracle PauliX application dependent on Fock Target ---
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
    print('The state of the ouput(in binary) is |' + str(q.split_register()) + '>')
    print('The target state(in binary) was |' + str(bin(target)[2:]) + '>')
    print('In Fock space this is |' + str(target) + '>')
    totalt = t.time() - init
    print('This took '+str(totalt)+' s to run')

main()
