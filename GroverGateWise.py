import numpy as np
import math as m
import time as t
import InOut as IO
from gatec import *
#from sparse import *
#from lazy import *

def numits(n):
    its = int((m.pi/4.0)*(2**n)**(1/2))

    return its

def findBinary(n, target):
    print('\nConverting Fock value to binary array...')
    ti = t.time()
    B = [int(x) for x in bin(target)[2:]]
    while len(B) != n:
        B = np.insert(B, 0, 0)
    Binaryform = B
    print('Binary array was formed in ' + str(t.time()-ti) + ' s')

    return Binaryform

def oracleX(n, Binaryform, x, I):
    print('\nAssigning PauliX gates to qubits for Oracle search...')
    ti = t.time()
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
    print('Assigning the PauliX gates took ' +str(t.time() - ti) + ' s')

    return Search

def grover(q, Search, cZ, H, X, its):
    #Create Superposition of states
    print('\nCreating superposition state...')
    t1 = t.time()
    q = H*q
    print('Creating superposition state took ' + str(t.time()-t1) + ' s')

    # # NOTE: Test Option - can comment in and use preformed O and D method
    # t2 = t.time()
    # print('\nForming the Oracle...')
    # Oracle = Search*cZ*Search
    # print('This took ' + str(t.time()-t2) + ' s')
    #
    # t3 = t.time()
    # print('\nForming the Diffusion matrix...')
    # Diffusion = H*X*cZ*X*H
    # print('This took ' + str(t.time()-t3) + ' s')

    #Grover's Iteration
    print('\nBeginning Grovers Iteration...')
    ti = t.time()
    for i in range(its):
        #Oracle
        q = Search*q
        q = cZ*q
        q = Search*q
        #Diffusion
        q = H*q
        q = X*q
        q = cZ*q
        q = X*q
        q = H*q

        # # NOTE: Test Option - commented in to use pre determined O and D
        # q = Oracle*q
        # q = Diffusion*q

        if i == 0:
            print('One Grover iteration took ' + str(t.time()-ti) + ' s')
    print('All of Grovers iteration took ' + str(t.time()-ti) + ' s')

    return q

def run(args, noise):
    # --- Reg size and target value ---
    n = args[0]
    target = args[1]

    # --- Timer Initialise ---
    t1 = t.time()

    # --- Initialised gates ---
    print('\nInitialising gates...')
    I = Identity()
    H = Hadamard(n)   #hadamard all
    X = PauliX(n)   #paulix all
    x = PauliX()
    z = PauliZ()
    cZ = Controlled(z, n)   #controlled z all
    #for noisy use
    if noise != 0:
        I = Noisy(I, noise)
        H = Noisy(H, noise)
        X = Noisy(X, noise)
        x = Noisy(x, noise)
        z = Noisy(z, noise)
        cZ = Noisy(cZ, noise)
    print('Gate initialisation took ' + str(t.time()-t1) + ' s')

    # --- Qreg formation ---
    print('\nForming quantum register...')
    t2 = t.time()
    q = Qubit(n)
    print('Quantum register formation took ' + str(t.time()-t2) + ' s')

    # --- Number of Iterations calculation ---
    its = numits(n)

    # --- Fock to Binary Array Conversion ---
    Binaryform = findBinary(n, target)

    # --- Oracle PauliX application dependent on Fock Target ---
    Search = oracleX(n, Binaryform, x, I)

    # --- Create Superposition and Grover's Iteration ---
    q = grover(q, Search, cZ, H, X, its)

    # --- Measure and Display ---
    q.measure()
    IO.printOut(q, target)
    print('\nThis took '+str(t.time()-t1)+' s to run\n')



def test(args):
    # --- Reg size and target value ---
    n = int(args[0])
    target = int(args[1])

    # --- Timer Initialise ---
    t1 = t.time()

    # --- Initialised gates ---
    I = Identity()
    H = Hadamard(n)   #hadamard all
    X = PauliX(n)   #paulix all
    x = PauliX()
    z = PauliZ()
    cZ = Controlled(z, n)   #controlled z all

    # --- Qreg formation ---
    q = Qubit(n)

    # --- Number of Iterations calculation ---
    its = numits(n)

    # --- Fock to Binary Array Conversion ---
    Binaryform = findBinary(n, target)

    # --- Oracle PauliX application dependent on Fock Target ---
    Search = oracleX(n, Binaryform, x, I)

    # --- Create Superposition and Grover's Iteration ---
    q = grover(q, Search, cZ, H, X, its)

    # --- Measure and Display ---
    q.measure()
    print('\nThis took '+str(t.time()-t1)+' s to run\n')

    return t.time()-t1
