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
            return Gate(self.array*other.array)

        elif (self.type=="Qubit") and (other.type=="Qubit"):
            assert self.array.shape==other.array.shape, "Qubit registers must have same size"
            return Gate(sp.outer(self.array,other.array))

        else:
            assert (self.type=="Gate") and (other.type=="Qubit"), "Gate must act on Qubit register"
            #assert self.array.shape[0]==other.array.shape[0], "Qubit register and gate must be of same size"
            return Qubit(other.array*self.array)



    def __str__(self):
        #String function for printing contents of SMatrix
        return str(self.array.toarray())


    def __len__(self):
        return int(np.log2(self.array.shape[1]))


    def __and__(self,other):
        #Magic method turning & into tensor products
        assert self.type==other.type, "Cannot tensor a Gate with a Qubit register"
        if (self.type=="Gate") and (other.type=="Gate"):
            
            return Gate(tensor_sparse_gate(self.array,other.array))
            #return Gate(sp.kron(self.array,other.array))
        elif (self.type=="Qubit") and (other.type=="Qubit"):
            return Qubit(tensor_sparse_qubit(self.array,other.array))
            #return Qubit(sp.kron(self.array,other.array))
    def ret(self):
        #returns array for plotting; not of type Qubit so works properly
        return self.array.toarray()

    def ret_mod(self):
        #Returns modulus squared of array i.e. Qubit register probabilities
        return np.abs(np.square(self.array.toarray()))[0]

class Hadamard(SMatrix):
    def __init__(self,n=1):
        SMatrix.__init__(self,"Gate")
        h = sp.bsr_matrix([[1,1],[1,-1]])
        hn = h
        for i in range(n-1):
            hn = tensor_sparse_gate(h,hn)
        hn = hn*(2**(-0.5*n))
        self.array = sp.bsr_matrix(hn)


class V(SMatrix):
    def __init__(self):
        SMatrix.__init__(self,"Gate")
        self.array = sp.bsr_matrix([[1,0],[0,1j]])

class Phase(SMatrix):
    def __init__(self,phase,n=1):
        SMatrix.__init__(self,"Gate")
        #self.phase = phase
        ph = sp.bsr_matrix([[1,0],[0,np.exp(1j*phase)]])
        phn = ph
        for i in range(n-1):
            ph = tensor_sparse_gate(ph,phn)
        self.array = sp.bsr_matrix(ph)

class Identity(SMatrix):
    def __init__(self,n=1):
        SMatrix.__init__(self,"Gate")
        self.array = sp.identity(2**n)

class PauliZ(SMatrix):
    def __init__(self):
        SMatrix.__init__(self,"Gate")
        self.array = sp.bsr_matrix([[1,0],
                                    [0,-1]])




class PauliX(SMatrix):
    def __init__(self,n=1):
        SMatrix.__init__(self,"Gate")
        self.array = sp.bsr_matrix(np.flipud(np.identity(2**n)))



##############################################################################################
#Control Gates

class Controlled(SMatrix):
    #General controlled gate. Takes any 1 qubit gate as input and makes a controlled version of that
    def __init__(self,other_gate,n=2):
        SMatrix.__init__(self,"Gate")
        self.array = sp.lil_matrix(sp.identity(2**n))
        t = other_gate.array.toarray()
        self.array[2**n-2,2**n-2] = t[0,0]
        self.array[2**n-1,2**n-1] = t[1,1]
        self.array[2**n-1,2**n-2] = t[1,0]
        self.array[2**n-2,2**n-1] = t[0,1]
        self.array = sp.bsr_matrix(self.array)




class CNot(SMatrix):
    def __init__(self,n=2):
        SMatrix.__init__(self,"Gate")
        self.array = sp.csr_matrix(sp.identity(2**n))
        self.array[2**n-2,2**n-2] = 0
        self.array[2**n-1,2**n-1] = 0
        self.array[2**n-1,2**n-2] = 1
        self.array[2**n-2,2**n-1] = 1
        self.array = sp.bsr_matrix(self.array)


class Toffoli(SMatrix):
    def __init__(self):
        SMatrix.__init__(self,"Gate")
        self.array = sp.bsr_matrix([[1,0,0,0,0,0,0,0],
                                    [0,1,0,0,0,0,0,0],
                                    [0,0,1,0,0,0,0,0],
                                    [0,0,0,1,0,0,0,0],
                                    [0,0,0,0,1,0,0,0],
                                    [0,0,0,0,0,1,0,0],
                                    [0,0,0,0,0,0,0,1],
                                    [0,0,0,0,0,0,1,0]])


class CPhase(SMatrix):
    def __init__(self,phase,n=2):
        SMatrix.__init__(self,"Gate")
        
        self.array = sp.csr_matrix(sp.identity(2**n),dtype=complex)
        self.array[2**n-1,2**n-1]=np.exp(1j*phase)
        self.array = sp.bsr_matrix(self.array,dtype=complex)
        #self.array = sp.bsr_matrix([[1,0,0,0],
        #                            [0,1,0,0],
        #                            [0,0,1,0],
        #                            [0,0,0,np.exp(1j*phase)]])


################################################################################
#Other useful gates

class Swap(SMatrix):
    def __init__(self,n=2,index1=0,index2=1):
        SMatrix.__init__(self,"Gate")

        self.array = sp.bsr_matrix(perm_matrix(n,index1,index2))

        

class Diffusion(SMatrix):
    def __init__(self,n):
        SMatrix.__init__(self,"Gate")
        N = 2**n
        c = 2.0/N
        self.array = sp.bsr_matrix(np.full((N, N), c) - np.identity(N))




class Oracle(SMatrix):
    def __init__(self,reg_size,target):
        SMatrix.__init__(self,"Gate")
        diags = np.ones(2**reg_size)
        #offsets = np.arange(0,2**reg_size,1)
        diags[target] = -1
        self.array = sp.csr_matrix(sp.identity(2**reg_size))
        self.array[target,target] = -1
        self.array = sp.bsr_matrix(self.array)
        #self.array = sp.dia_matrix(diags)

class Gate(SMatrix):
    #Generic gate class - used as output for multiplication or tensor of other gates
    def __init__(self,data):
        SMatrix.__init__(self,"Gate")
        #assert (len(data[0])&(len(data[0])-1)==0) and (len(data[1])&(len(data[1])-1)==0) and (len(data[0])==len(data[1])),"Gate must be square matrix of size 2**n"
        self.array = sp.bsr_matrix(data)

class Noisy(SMatrix):
    def __init__(self,other_gate,a=0.1):
        #generate a noisy version a gate
        SMatrix.__init__(self,"Gate")
        self.array = other_gate.array
        self.array = self.array + sp.random(self.array.shape[0],self.array.shape[0],density=a)
        self.array = self.array + 1j*sp.random(self.array.shape[0],self.array.shape[0],density=a)


class Qubit(SMatrix):
    #Class for Qubit
    def __init__(self,data,fock=0):
        SMatrix.__init__(self,"Qubit")
        if type(data) is int:
            self.array = np.zeros(2**data)
            self.array[fock] = 1
            self.array = sp.bsr_matrix(self.array) 
        else:
            self.array = sp.bsr_matrix(data)
        


        #print(self.array.shape)
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
        data = self.array.toarray()[0]
        pos = np.arange(len(data))
        #print(pos)
        probs = np.abs(np.square(data))
        #If probs is not normalised (usually due to rounding errors), re-normalise
        probs = probs/np.sum(probs)
        #print(probs)
        dist = stats.rv_discrete(values=(pos,probs))
        self.array = np.zeros(data.shape)
        self.array[dist.rvs()] = 1
        self.array = sp.bsr_matrix(self.array)

    def measure_cheat(self):
        #Measure but ignore 0 state, for debugging shors
        data = self.array.toarray()[0]
        pos = np.arange(len(data))
        #print(pos)
        probs = np.abs(np.square(data))
        probs[0] = 0
        #If probs is not normalised (usually due to rounding errors), re-normalise
        probs = probs/np.sum(probs)
        #print(probs)
        dist = stats.rv_discrete(values=(pos,probs))
        self.array = np.zeros(data.shape)
        self.array[dist.rvs()] = 1
        self.array = sp.bsr_matrix(self.array)

    def split_register(self):
        #Only run after measured. returns individual qubit values
        t = self.array.toarray()[0]
        outs = np.arange(0,len(t),1)
        res = np.array(np.sum(outs*t.astype(int)))
        return np.binary_repr(res)




#def main():
    #h = Hadamard(10)
    #q = Qubit([1,0])
    #v = V()
    #i = Identity()
    #print(v)
    #print(h)
    #print(q)
    #print(i&h)
    #print(v*h*q)
    #q1 = h*q
    #print(q1)
    #print(v*q1)
    #g = v*h
    #print(g*q)

    #print(Hadamard(2))

#main()