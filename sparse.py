import numpy as np
from scipy import stats
from scipy import sparse as sp
import math
import cmath
from utilities import *


#Sparse matrix optimisation of low level bits - AR
#use either csr or bsr matrix



class SMatrix:
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
            return Gate(sp.matmul(self.array,other.array))

        elif (self.type=="Qubit") and (other.type=="Qubit"):
            assert self.array.shape==other.array.shape, "Qubit registers must have same size"
            return Gate(sp.outer(self.array,other.array))

        else:
            assert (self.type=="Gate") and (other.type=="Qubit"), "Gate must act on Qubit register"
            assert self.array.shape[0]==other.array.shape[0], "Qubit register and gate must be of same size"
            return Qubit(sp.matmul(self.array,other.array))



    def __str__(self):
        #String function for printing contents of SMatrix
        return str(self.array.toarray())


    def __len__(self):
        return int(np.log2(len(self.array)))


    def __and__(self,other):
        #Magic method turning & into tensor products
        assert self.type==other.type, "Cannot tensor a Gate with a Qubit register"
        if (self.type=="Gate") and (other.type=="Gate"):
            
            #return Gate(tensor(self.array,other.array))
            return Gate(sp.kron(self.array,other.array))
        elif (self.type=="Qubit") and (other.type=="Qubit"):
            #return Qubit(tensor(self.array,other.array))
            return Qubit(sp.kron(self.array,other.array))
    def ret(self):
        #returns array for plotting; not of type Qubit so works properly
        return self.array


class Hadamard(SMatrix):
    def __init__(self,n=1):
        SMatrix.__init__(self,"Gate")
        h = np.array([[1,1],[1,-1]])
        hn = h
        for i in range(n-1):
            hn = tensor(h,hn)
        hn = hn*(2**(-0.5*n))
        self.array = sp.bsr_matrix(np.array(hn))

class Diffusion(SMatrix):
    def __init__(self,n):
        SMatrix.__init__(self,"Gate")
        N = 2**n
        c = 2.0/N
        self.array = sp.bsr_matrix(np.full((N, N), c) - np.identity(N))

class V(SMatrix):
    def __init__(self):
        SMatrix.__init__(self,"Gate")
        self.array = sp.bsr_matrix(np.array([[1,0],[0,1j]]))

class Phase(SMatrix):
    def __init__(self,phase,n=1):
        SMatrix.__init__(self,"Gate")
        self.phase = phase
        ph = np.array([[1,0],[0,np.exp(1j*phase)]])
        phn = ph
        for i in range(n-1):
            ph = tensor(ph,phn)
        self.array = sp.bsr_matrix(ph)

class Identity(SMatrix):
    def __init__(self,n=1):
        SMatrix.__init__(self,"Gate")
        self.array = sp.bsr_matrix(np.identity(2**n))


class PauliX(SMatrix):
    def __init__(self,n=1):
        SMatrix.__init__(self,"Gate")
        self.array = sp.bsr_matrix(np.flipud(np.identity(2**n)))
# 2 Qubit Gates

class CNot(SMatrix):
    def __init__(self):
        SMatrix.__init__(self,"Gate")
        self.array = sp.bsr_matrix(np.array([[1,0,0,0],
                               [0,1,0,0],
                               [0,0,0,1],
                               [0,0,1,0]]))

class CPhase(SMatrix):
    def __init__(self,phase):
        SMatrix.__init__(self,"Gate")
        self.array = sp.bsr_matrix(np.array([[1,0,0,0],
                               [0,1,0,0],
                               [0,0,1,0],
                               [0,0,0,np.exp(1j*phase)]]))

class Swap(SMatrix):
    def __init__(self,n=2,index1=0,index2=1):
        SMatrix.__init__(self,"Gate")

        self.array = sp.bsr_matrix(perm_matrix(n,index1,index2))

        #self.array = np.array([[1,0,0,0],
        #                       [0,0,1,0],
        #                       [0,1,0,0],
        #                       [0,0,0,1]])

# 3 qubit gates

class Toffoli(SMatrix):
    def __init__(self):
        SMatrix.__init__(self,"Gate")
        self.array = sp.bsr_matrix(np.array([[1,0,0,0,0,0,0,0],
                               [0,1,0,0,0,0,0,0],
                               [0,0,1,0,0,0,0,0],
                               [0,0,0,1,0,0,0,0],
                               [0,0,0,0,1,0,0,0],
                               [0,0,0,0,0,1,0,0],
                               [0,0,0,0,0,0,0,1],
                               [0,0,0,0,0,0,1,0]]))




class Oracle(SMatrix):
    def __init__(self,reg_size,target):
        SMatrix.__init__(self,"Gate")
        self.array = sp.bsr_matrix(np.identity(2**reg_size))
        self.array[target,target] = -1

class Gate(SMatrix):
    #Generic gate class - used as output for multiplication or tensor of other gates
    def __init__(self,data):
        SMatrix.__init__(self,"Gate")
        #assert (len(data[0])&(len(data[0])-1)==0) and (len(data[1])&(len(data[1])-1)==0) and (len(data[0])==len(data[1])),"Gate must be square matrix of size 2**n"
        self.array = sp.bsr_matrix(data)



class Qubit(SMatrix):
    #Class for Qubit
    def __init__(self,data):
        SMatrix.__init__(self,"Qubit")
        assert (len(data)&(len(data)-1)==0),"Qubit register length must be a power of 2"
        self.array = sp.bsr_matrix(data)
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
        self.array = sp.bsr_matrix(np.divide(self.array,a))

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

        outs = np.arange(1,len(self.array)+1,1)
        res = np.array(np.sum(outs*self.array.astype(int)))
        return np.binary_repr(res)




def main():
    h = Hadamard(1)
    print(h&h)
    print(Hadamard(2))

main()