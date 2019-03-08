import numpy as np
import math as m
import time as t
#import InOut as IO
#from gatec import *
#from sparse import *
from lazy import *
#from lazzz import *

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
    #print(H)
    #print(q)
    q = H*q
    print('Creating superposition state took ' + str(t.time()-t1) + ' s')
    #Grover's Iteration
    print('\nBeginning Grovers Iteration...')
    ti = t.time()
    for i in range(its):
        q = (Search&cZ&Search&H&X&cZ&X&H)*q
       # print((q))
        #print("Search")
        #q = Search*q
        #THIS IS THE PART WHERE PRINTING /EVAL STOPS WORKING
        #print(type(q))
        #print(q)
       # print("cz")
        #q = cZ*q
        #print("Search")
        #q = Search*q
        #print("H")
        #q = H*q
        #print("X")
        #q = X*q
        #print("cz")


        #q = cZ*q
        #q = X*q
        #q = H*q
        #print(q)
        if i == 0:
            print('One Grover iteration took ' + str(t.time()-ti) + ' s')
    print('All of Grovers iteration took ' + str(t.time()-ti) + ' s')

    return q

def main():
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
    #print(q)
    # --- Number of Iterations calculation ---
    its = int((m.pi/4.0)*(2**n)**(1/2))

    # --- Fock to Binary Array Conversion ---
    Binaryform = findBinary(n, target)
    #print(q)
    # --- Oracle PauliX application dependent on Fock Target ---
    Search = oracleX(n, Binaryform, x, I)
    #print(q)
    # --- Create Superposition and Grover's Iteration ---
    q = grover(q, Search, cZ, H, X, its)
    # --- Measure and Display ---
    q.measure()
    print('\nThe state of the ouput(in binary) is |' + str(q.split_register()) + '>')
    print('The target state(in binary) was |' + str(bin(target)[2:]) + '>')
    print('In Fock space this is |' + str(target) + '>')
    print('\nThis took '+str(t.time()-t1)+' s to run\n')

main()
