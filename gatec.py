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

        #multiplication order:
        #    Gate*Qubit
        #    Gate*Scalar
        #    Qubit*Scalar
        #
        if (type(other)==int):
            if (self.type=="Gate"):
                return Gate(self.array*other)
            else:
                return Qubit(self.array*other)

        elif (self.type=="Gate") and (other.type=="Gate"):
            assert self.array.shape==other.array.shape, "Gates must be of same size"
            return Gate(np.matmul(self.array,other.array))

        elif (self.type=="Qubit") and (other.type=="Qubit"):
            assert self.array.shape==other.array.shape, "Qubit registers must have same size"
            return Gate(np.outer(self.array,other.array))
        
        else:
            assert (self.type=="Gate") and (other.type=="Qubit"), "Gate must act on Qubit register"
            assert self.array.shape[0]==other.array.shape[0], "Qubit register and gate must be of same size"
            return Qubit(np.matmul(self.array,other.array))

        

    def __str__(self):
        #String function for printing contents of QMatrix
        return str(self.array)


    def get_size(self):
        return np.log2(len(self.array))


    def __and__(self,other):
        #Magic method turning & into tensor products
        assert self.type==other.type, "Cannot tensor a Gate with a Qubit register"
        if (self.type=="Gate") and (other.type=="Gate"):
            return Gate(tensor(self.array,other.array))
        elif (self.type=="Qubit") and (other.type=="Qubit"):
            return Qubit(tensor(self.array,other.array))



#------------------------------------------------------
#---some gates and qubits under QMatrix parent class.
#------------------------------------------------------

# Single Qubit Gates

class Hadamard(QMatrix):
    def __init__(self,n=1):
        QMatrix.__init__(self,"Gate")
        h = np.array([[1,1],[1,-1]])
        hn = h
        for i in range(n-1):
            hn = tensor(h,hn)
        hn = hn*(2**(-0.5*n))
        self.array = np.array(hn)


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
        assert (len(data[0])&(len(data[0])-1)==0) and (len(data[1])&(len(data[1])-1)==0) and (len(data[0])==len(data[1])),"Gate must be square matrix of size 2**n"
        self.array = np.array(data)



class Qubit(QMatrix):
    #Class for Qubit
    def __init__(self,data):
        QMatrix.__init__(self,"Qubit")
        assert (len(data)&(len(data)-1)==0),"Qubit register length must be a power of 2"
        self.array = np.array(data)


    



def main():

    #Just about manages 10 qubits, but will get uncomfortably slow beyond that

    q = Qubit([1,0])
    
    

    #h4 = Hadamard(4)
    h2 = Hadamard(2)
    q2 = h2*(q&q)



    #h = q&h2
    
    
    print(q2)
    print(h2*2)

    #print(f1)
    #print(f2)
    #b = V()
    #c = a*b
    #print(c.array)

main()