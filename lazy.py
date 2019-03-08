import numpy as np
from scipy import stats
#from scipy import sparse as sp
import lazyarray
from lazyarray import larray
from lazyarray import matmul
import math
import cmath
from utilities import *


#Sparse matrix optimisation of low level bits - AR
#use either csr or bsr matrix



class LMatrix:
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
            return Gate(lazy_mul(self.array,other.array))

        elif (self.type=="Qubit") and (other.type=="Qubit"):
            assert self.array.shape==other.array.shape, "Qubit registers must have same size"
            return Gate(lazy_mul(self.array,other.array))

        else:
            assert (self.type=="Gate") and (other.type=="Qubit"), "Gate must act on Qubit register"
            #assert self.array.shape[0]==other.array.shape[0], "Qubit register and gate must be of same size"

            return Qubit(lazy_qub_mul(self.array,other.array))
        
       


    def __str__(self):
        #String function for printing contents of LMatrix
        return str(self.array.evaluate())


    def __len__(self):
        return int(np.log2(self.array.shape[1]))


    def __and__(self,other):
        #Magic method turning & into tensor products
        assert self.type==other.type, "Cannot tensor a Gate with a Qubit register"
        if (self.type=="Gate") and (other.type=="Gate"):
            
            return Gate(tensor_lazy(self.array,other.array))
            #return Gate(sp.kron(self.array,other.array))
        elif (self.type=="Qubit") and (other.type=="Qubit"):
            return Qubit(tensor_lazy(self.array,other.array))
            #return Qubit(sp.kron(self.array,other.array))
    def ret(self):
        #returns array for plotting; not of type Qubit so works properly
        return self.array.toarray()


class Hadamard(LMatrix):
    def __init__(self,n=1):
        LMatrix.__init__(self,"Gate")
        h = larray([[1,1],[1,-1]])
        hn = h
        for i in range(n-1):
            hn = tensor_lazy(h,hn)
        hn = hn*(2**(-0.5*n))
        self.array = larray(hn)


class Diffusion(LMatrix):
    def __init__(self,n):
        LMatrix.__init__(self,"Gate")
        N = 2**n
        c = 2.0/N
        self.array = larray(np.full((N, N), c) - np.identity(N))

class V(LMatrix):
    def __init__(self):
        LMatrix.__init__(self,"Gate")
        self.array = larray([[1,0],[0,1j]])

class Phase(LMatrix):
    def __init__(self,phase,n=1):
        LMatrix.__init__(self,"Gate")
        self.phase = phase
        ph = larray([[1,0],[0,np.exp(1j*phase)]])
        phn = ph
        for i in range(n-1):
            ph = tensor_lazy(ph,phn)
        self.array = larray(ph)

class Identity(LMatrix):
    def __init__(self,n=1):
        LMatrix.__init__(self,"Gate")
        self.array = larray(np.identity(2**n))


class PauliX(LMatrix):
    def __init__(self,n=1):
        LMatrix.__init__(self,"Gate")
        self.array = larray(np.flipud(np.identity(2**n)))
        
class PauliZ(LMatrix):
    def __init__(self):
        LMatrix.__init__(self,"Gate")
        self.array = larray([[1,0],
                             [0,-1]])
# 2 Qubit Gates



class CNot(LMatrix):
    def __init__(self,n=2):
        SMatrix.__init__(self,"Gate")
        self.array = sp.csr_matrix(sp.identity(2**n))
        self.array[2**n-2,2**n-2] = 0
        self.array[2**n-1,2**n-1] = 0
        self.array[2**n-1,2**n-2] = 1
        self.array[2**n-2,2**n-1] = 1
        self.array = larray(self.array)
        #LMatrix.__init__(self,"Gate")
        #self.array = larray([[1,0,0,0],
        #                            [0,1,0,0],
        #                            [0,0,0,1],
        #                            [0,0,1,0]])

class CPhase(LMatrix):
    def __init__(self,phase):
        LMatrix.__init__(self,"Gate")
        self.array = larray([[1,0,0,0],
                                    [0,1,0,0],
                                    [0,0,1,0],
                                    [0,0,0,np.exp(1j*phase)]])

class Swap(LMatrix):
    def __init__(self,n=2,index1=0,index2=1):
        LMatrix.__init__(self,"Gate")

        self.array = larray(perm_matrix(n,index1,index2))

        

# 3 qubit gates

class Toffoli(LMatrix):
    def __init__(self):
        LMatrix.__init__(self,"Gate")
        self.array = larray([[1,0,0,0,0,0,0,0],
                                    [0,1,0,0,0,0,0,0],
                                    [0,0,1,0,0,0,0,0],
                                    [0,0,0,1,0,0,0,0],
                                    [0,0,0,0,1,0,0,0],
                                    [0,0,0,0,0,1,0,0],
                                    [0,0,0,0,0,0,0,1],
                                    [0,0,0,0,0,0,1,0]])



class Oracle(LMatrix):
    def __init__(self,reg_size,target):
        LMatrix.__init__(self,"Gate")
        diags = np.ones(2**reg_size)
        #offsets = np.arange(0,2**reg_size,1)
        diags[target] = -1
        self.array = sp.csr_matrix(sp.identity(2**reg_size))
        self.array[target,target] = -1
        self.array = larray(self.array)
        #self.array = sp.dia_matrix(diags)

class Gate(LMatrix):
    #Generic gate class - used as output for multiplication or tensor of other gates
    def __init__(self,data):
        LMatrix.__init__(self,"Gate")
        #assert (len(data[0])&(len(data[0])-1)==0) and (len(data[1])&(len(data[1])-1)==0) and (len(data[0])==len(data[1])),"Gate must be square matrix of size 2**n"
        self.array = data

class Controlled(LMatrix):
    #General controlled gate. Takes any 1 qubit gate as input and makes a controlled version of that
    def __init__(self,other_gate,n=2):
        LMatrix.__init__(self,"Gate")
        self.array = (np.identity(2**n))
        t = other_gate.array.evaluate()
        self.array[2**n-2,2**n-2] = t[0,0]
        self.array[2**n-1,2**n-1] = t[1,1]
        self.array[2**n-1,2**n-2] = t[1,0]
        self.array[2**n-2,2**n-1] = t[0,1]
        self.array = larray(self.array)
        #print(self.array.evaluate())

class Qubit(LMatrix):
    #Class for Qubit
    def __init__(self,data,fock=0):
        LMatrix.__init__(self,"Qubit")
        if type(data) is int:
            self.array = np.zeros(2**data)
            self.array[fock] = 1
            self.array = [self.array]
            self.array = larray(self.array)
        else:
            self.array = data


    def normalise(self):
        div = np.sqrt(np.sum(np.square(self.array)))
        a = np.empty(len(self.array))
        a.fill(div)
        self.array = larray(np.divide(self.array,a))

    def measure(self):
        #method to collapse qubit register into 1 state.
        print (self.array.evaluate())
        data = self.array.evaluate()
        print(type(data))
        pos = np.arange(len(data))
        #print(pos)
        probs = np.abs(np.square(data))
        #If probs is not normalised (usually due to rounding errors), re-normalise
        probs = probs/np.sum(probs)
        #print(probs)
        dist = stats.rv_discrete(values=(pos,probs))
        self.array = np.zeros(data.shape)
        self.array[dist.rvs()] = 1
        return self.array

    def split_register(self):
        #Only run after measured. returns individual qubit values

        outs = np.arange(0,len(self.array),1)
        res = np.array(np.sum(outs*self.array.astype(int)))
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
