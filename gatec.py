import numpy as np
from scipy import stats
import math
import cmath
from utilities import *






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


    def __len__(self):
        return int(np.log2(len(self.array)))


    def __and__(self,other):
        #Magic method turning & into tensor products
        assert self.type==other.type, "Cannot tensor a Gate with a Qubit register"
        if (self.type=="Gate") and (other.type=="Gate"):
            return Gate(tensor(self.array,other.array))
        elif (self.type=="Qubit") and (other.type=="Qubit"):
            return Qubit(tensor(self.array,other.array))

    def ret(self):
        #returns array for plotting; not of type Qubit so works properly
        return self.array

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

class Diffusion(QMatrix):
    def __init__(self,n):
        QMatrix.__init__(self,"Gate")
        N = 2**n
        c = 2.0/N
        self.array = np.full((N, N), c) - np.identity(N)

class V(QMatrix):
    def __init__(self):
        QMatrix.__init__(self,"Gate")
        self.array = np.array([[1,0],[0,1j]])

class Phase(QMatrix):
    def __init__(self,phase,n=1):
        QMatrix.__init__(self,"Gate")
        self.phase = phase
        ph = np.array([[1,0],[0,np.exp(1j*phase)]])
        phn = ph
        for i in range(n-1):
            ph = tensor(ph,phn)
        self.array = ph

class Identity(QMatrix):
    def __init__(self,n=1):
        QMatrix.__init__(self,"Gate")
        self.array = np.identity(2**n)

class PauliZ(QMatrix):
    def __init__(self):
        QMatrix.__init__(self,"Gate")
        self.array = np.array([[1,0],[0,-1]])

class PauliX(QMatrix):
    def __init__(self,n=1):
        QMatrix.__init__(self,"Gate")
        self.array = np.flipud(np.identity(2**n))


##############################################################################################
#Control Gates


class Controlled(QMatrix):
    #General controlled gate. Takes any 1 qubit gate as input and makes a controlled version of that
    def __init__(self,other_gate,n=2):
        QMatrix.__init__(self,"Gate")
        self.array = np.identity(2**n)
        self.array[2**n-2,2**n-2] = other_gate.array[0,0]
        self.array[2**n-1,2**n-1] = other_gate.array[1,1]
        self.array[2**n-1,2**n-2] = other_gate.array[1,0]
        self.array[2**n-2,2**n-1] = other_gate.array[0,1]
        


class CNot(QMatrix):
    def __init__(self,n=2):
        QMatrix.__init__(self,"Gate")
        self.array = np.identity(2**n)
        self.array[2**n-2,2**n-2] = 0
        self.array[2**n-1,2**n-1] = 0
        self.array[2**n-1,2**n-2] = 1
        self.array[2**n-2,2**n-1] = 1
        
class Toffoli(QMatrix):
    def __init__(self):
        QMatrix.__init__(self,"Gate")
        self.array = np.array([[1,0,0,0,0,0,0,0],
                               [0,1,0,0,0,0,0,0],
                               [0,0,1,0,0,0,0,0],
                               [0,0,0,1,0,0,0,0],
                               [0,0,0,0,1,0,0,0],
                               [0,0,0,0,0,1,0,0],
                               [0,0,0,0,0,0,0,1],
                               [0,0,0,0,0,0,1,0]])


class CPhase(QMatrix):
    def __init__(self,phase,n=2):
        QMatrix.__init__(self,"Gate")
        self.array = np.identity(2**n)
        self.array[2**n-1,2**n-1]=np.exp(1j*phase)

        #self.array = np.array([[1,0,0,0],
        #                       [0,1,0,0],
        #                       [0,0,1,0],
        #                       [0,0,0,np.exp(1j*phase)]])

################################################################
#Other useful gates



class Swap(QMatrix):
    def __init__(self,n=2,index1=0,index2=1):
        QMatrix.__init__(self,"Gate")

        self.array = perm_matrix(n,index1,index2)

        #self.array = np.array([[1,0,0,0],
        #                       [0,0,1,0],
        #                       [0,1,0,0],
        #                       [0,0,0,1]])

# 3 qubit gates





class Oracle(QMatrix):
    def __init__(self,reg_size,target):
        QMatrix.__init__(self,"Gate")
        self.array = np.identity(2**reg_size)
        self.array[target,target] = -1

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
        #catches and normalises unnormalised qubits. good for testing but shouldnt be needed in the end
        #Causes errors, commented out for now
        #if 0.9999999 < (np.sum(np.square(self.array))) < 1.00000001:
        #    pass
        #else:
        #    self.normalise()

    def normalise(self):
        div = np.sqrt(np.sum(np.square(self.array)))
        a = np.empty(len(self.array))
        a.fill(div)
        self.array = np.divide(self.array,a)

    def measure(self):
        #method to collapse qubit register into 1 state.
        pos = np.arange(len(self.array))
        probs = np.abs(np.square(self.array))
        #If probs is not normalised (usually due to rounding errors), re-normalise
        probs = probs/np.sum(probs)
        dist = stats.rv_discrete(values=(pos,probs))
        self.array = np.zeros(self.array.shape)
        self.array[dist.rvs()] = 1
        return self.array

    def split_register(self):
        #Only run after measured. returns individual qubit values

        outs = np.arange(0,len(self.array),1)
        res = np.array(np.sum(outs*self.array.astype(int)))
        return np.binary_repr(res)
