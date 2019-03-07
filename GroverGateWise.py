import numpy as np
import math as m
import time as t
import InOut as IO
#from gatec import *
from sparse import *

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

    #Grover's Iteration
    print('\nBeginning Grovers Iteration...')
    ti = t.time()
    for i in range(its):
        q = Search*q
        q = cZ*q
        q = Search*q

        q = H*q
        q = X*q
        q = cZ*q
        q = X*q
        q = H*q
        if i == 0:
            print('One Grover iteration took ' + str(t.time()-ti) + ' s')
    print('All of Grovers iteration took ' + str(t.time()-ti) + ' s')

    return q

def main():
    # --- How to run ---
    io = IO.start()
    n = io[0]
    target = io[1]

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
    print('Gate initialisation took ' + str(t.time()-t1) + ' s')

    # --- Qreg formation ---
    print('\nForming quantum register...')
    t2 = t.time()
    q = Qubit(n)
    print('Quantum register formation took ' + str(t.time()-t2) + ' s')

    # --- Number of Iterations calculation ---
    its = int((m.pi/4.0)*(2**n)**(1/2))

    # --- Fock to Binary Array Conversion ---
    Binaryform = findBinary(n, target)

    # --- Oracle PauliX application dependent on Fock Target ---
    Search = oracleX(n, Binaryform, x, I)

    # --- Create Superposition and Grover's Iteration ---
    q = grover(q, Search, cZ, H, X, its)

    # --- Measure and Display ---
    q.measure()
    print('\nThe state of the ouput(in binary) is |' + str(q.split_register()) + '>')
    print('The target state(in binary) was |' + str(bin(target)[2:]) + '>')
    print('In Fock space this is |' + str(target) + '>')
    print('\nThis took '+str(t.time()-t1)+' s to run\n')

main()
