import numpy as np
import math as m
import time as t
import InOut as IO
#from gatec import *
from sparse import *

def main():
    t1 = t.time()
    print("\nInitialise things:")
    # --- Number of qubits and Target Fock value ---
    n = 4
    target = 15
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
    t2 = t.time()

    print(t2-t1)
    print("\nSome initial calculations")
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
    t3 = t.time()
    print(t3-t2)
    print("\nBinary looping")
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
    t4 = t.time()
    print(t4-t3)
    # print("\nCalculate oracle and diffusion gates")
    # # --- Oracle and Diffusion Gate Calculation ---
    # Oracle = Search*cZ*Search
    # tcheck1 = t.time()
    # tcheck = tcheck1 - t4
    # print('Time for Orcacle calc '+str(tcheck))
    # Diffusion = H*X*cZ*X*H
    # tcheck2 = t.time()
    # tcheck = tcheck2 - tcheck1
    # print('Time for Diffusion calc '+str(tcheck))

    IO.Display(Search)
    IO.Display(cZ)
    # --- Initialising Superposition State ---
    q = H*q

    # --- Grover's Iteration ---
    t5 = t.time()
    print(t5-t4)
    print("\nRun grover's")
    for i in range(its):
        tloop = t.time()
        q = Search*q
        q = cZ*q
        q = Search*q

        q = H*q
        q = X*q
        q = cZ*q
        q = X*q
        q = H*q
        print('1 GI time: ' + str(t.time()-tloop))

    # --- Measure and Display ---
    q.measure()
    t6 = t.time()
    print(t6-t5)

<<<<<<< HEAD





=======
>>>>>>> 5db1831f69580965df4f491fdbce321901376fac
    print('\nThe binary value of the ouput is ' + str(q.split_register()))
    #totalt = t.time() - init
    print("\nTotal time:")
    print(t6-t1)
    #print('This took '+str(totalt)+' s to run')

main()
