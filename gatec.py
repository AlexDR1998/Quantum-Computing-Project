import numpy as np
import math
import cmath

#qubit class?
class Qubit:
    def __init__(self,state0):
        self.state = np.array([[state0[0]],[state0[1]]])


class Gate:
#parent class for all gates

    def __init__(self, rows, columns):
        self.type = "Gate" #test, delete later

    def __mul__(self,other):
        #magic method turning all * into matrix multiplication
        return np.matmul(self.array,other)

#various gates under Gate parent class.

class Hadamard(Gate):
    def __init__(self):
        Gate.__init__(self,0,0)
        self.array = np.array([[1,1],[1,-1]])*(2**(-0.5))

class ControlledV(Gate):
    def __init__(self):
        Gate.__init__(self,0,0)
        self.array = np.array([[1,0],[0,1j]])

class Phase_shift(Gate):
    def __init__(self,phase):
        Gate.__init__(self,0,0)
        self.phase = phase
        self.array = np.array([[1,0],[0,cmath.exp(1j*phase)]])
