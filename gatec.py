import numpy as np
import math
import cmath
from utilities import tensor




class QMatrix:
# abstract parent class for all gates and qubits

    def __init__(self, typ):
        self.type = str(typ) #Differentiates between Qubits and Gates

    def __mul__(self,other):
        #magic method turning all * into matrix multiplication
        

        if (self.type=="Gate") and (other.type=="Gate"):
            return Gate(np.matmul(self.array,other.array))
        else:
            return Qubit(np.matmul(self.array,other.array))

        
    def __str__(self):
        #String function for printing contents of QMatrix
        return str(self.array)

    def __and__(self,other):
        #Magic method turning & into tensor products

        if (self.type=="Gate") and (other.type=="Gate"):
            return Gate(tensor(self.array,other.array))
        elif (self.type=="Qubit") and (other.type=="Qubit"):
            return Qubit(tensor(self.array,other.array))



#------------------------------------------------------
#---some gates and qubits under QMatrix parent class.
#------------------------------------------------------

# Single Qubit Gates

class Hadamard(QMatrix):
    def __init__(self):
        QMatrix.__init__(self,"Gate")
        self.array = np.array([[1,1],[1,-1]])*(2**(-0.5))

class V(QMatrix):
    def __init__(self):
        QMatrix.__init__(self,"Gate")
        self.array = np.array([[1,0],[0,1j]])

class Phase(QMatrix):
    def __init__(self,phase):
        QMatrix.__init__(self,"Gate")
        self.phase = phase
        self.array = np.array([[1,0],[0,np.exp(1j*phase)]])




# 2 Qubit Gates

class CNot(QMatrix):
    def __init__(self):
        QMatrix.__init__(self,"Gate")
        self.array = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])

class CPhase(QMatrix):
    def __init__(self,phase):
        QMatrix.__init__(self,"Gate")
        self.array = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,np.exp(1j*phase)]])


class Gate(QMatrix):
    #Generic gate class - used as output for multiplication or tensor of other gates
    def __init__(self,data):
        QMatrix.__init__(self,"Gate")
        self.array = np.array(data)



class Qubit(QMatrix):
    #Class for Qubit
    def __init__(self,data):
        QMatrix.__init__(self,"Qubit")
        assert (len(data) & (len(data) - 1)==0),"Qubit register length must be a power of 2"
        self.array = np.array(data)


    def get_size(self):
        return np.log2(len(self.array))



def main():



    #q1 = a*q
    q1 = Qubit([1/np.sqrt(2),1/np.sqrt(2)])
    q2 = Qubit([0,1])

    q12 = q1&q2


    h1 = Hadamard()
    h2 = Hadamard()

    h12 = h1&h2


    f1 = q12*h12

    f2 = (q1*h1)&(q2*h2)

    print(f1)
    print(f2)
    #b = V()
    #c = a*b
    #print(c.array)

main()