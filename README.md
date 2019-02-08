# Quantum-Computing-Project

ALright kids this is your handy guide to using our low level quantum stuff. relevant 7/2/2019

1. What you need:
to use this stuff you need to from gatec.py import * into your code
its also a good idea to get numpy just in case

2. Creating Objects:
before you calculate your algorithm you need to create the gates and qubits you are going to use. They are all matrices so share some similar properties. Calling print() on any of these will return the array, and len() will return the number of qubits it operates on (not the actual size of the matrix)

qubits:
q = Qubit([1d array of values]) 
q.measure() measures the qubit


gates (important ones: have a look for others if you need)
Hadamard: H = Hadamard()
Identity: I = Identity(n) where n is number of qubits you are operating on, default is 1
CPhase: B = CPhase(phase)
Pauli X: px = PauliX()

You can also make your own gate (probably useful for oracle?) with:
mygate = Gate([2d array of values])

3. Operations:
tensor product: & symbol between matrices
matrix multiplication: * 
assert errors should be helpful enough here, remember laws of matrix multiplication

remember matrices are applied from R to L (according to J)

if this very sparse guide is not helpful,harrass me on the groupchat.

Gralster xoxox




