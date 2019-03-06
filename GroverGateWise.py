import numpy as np
import math as m
import time as t
import InOut as IO
#from gatec import *
from sparse import *

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
    t1 = t.time()

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
    print('\nForming quantum register...')
    t2 = t.time()
    q = q0&q0
    if n > 2:
        for i in range(n-2):
            q = q0&q
    print('Quantum register formation took ' + str(t.time()-t2) + ' s')

    # --- Number of Iterations calculation ---
    its = int((m.pi/4.0)*(N)**(1/2))

    # --- Fock to Binary Array Conversion ---
    print('\nConverting Fock value to binary array...')
    t3 = t.time()
    B = [int(x) for x in bin(target)[2:]]
    while len(B) != n:
        B = np.insert(B, 0, 0)
    Binaryform = B
    print('Binary array was formed in ' + str(t.time()-t3) + ' s')

    # --- Oracle PauliX application dependent on Fock Target ---
    print('\nAssigning PauliX gates to qubits for Oracle search...')
    t4 = t.time()
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
    print('Assigning the PauliX gates took ' +str(t.time() - t4) + ' s')

    # --- Initialising Superposition State ---
    print('\nCreating superposition state...')
    t5 = t.time()
    q = H*q
    print('Creating superposition state took ' + str(t.time()-t5) + ' s')

    # --- Grover's Iteration ---
    print('\nBeginning Grovers Iteration...')
    t6 = t.time()
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
            print('One Grover iteration took ' + str(t.time()-t6) + ' s')
    print('All of Grovers iteration took ' + str(t.time()-t6) + ' s')

    # --- Measure and Display ---
    q.measure()
    print('\nThe state of the ouput(in binary) is |' + str(q.split_register()) + '>')
    print('The target state(in binary) was |' + str(bin(target)[2:]) + '>')
    print('In Fock space this is |' + str(target) + '>')
    print('\nThis took '+str(t.time()-t1)+' s to run\n')

main()
