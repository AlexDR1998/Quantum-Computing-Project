import numpy as np
from scipy import stats
import math
import cmath
from utilities import *

"""
Dense implementation of low level (gates and Qubits) as matrices
"""

class QMatrix:
    """Abstract parent class for all quantum objects
    """

    def __init__(self, typ):
        """Set the type of object: Gate or Qubit
        """
        self.type = str(typ)

    def __mul__(self,other):
        """Method to perform matrix multiplication. Combines gates sequentially
        and acts gates on registers. Multiplied objects must have the same size

        """
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

    def get(self,n):
        return str(self.array[n])

    def __str__(self):
        """Returns string of contents of quantum object. Not 'realistic', as
        quantum states cannot be fully observed without collapsing qubit register
        to 1 state. For that, use Qubit.measure()
        """
        return str(self.array)


    def __len__(self):
        """Returns the size of gate or qubit register, in number of qubits
        """
        return int(np.log2(len(self.array)))


    def __and__(self,other):
        """Method to perform tensor products of quantum objects. Combines
        gates in parallel. Combines qubits into qubit registers, or qubit registers
        into bigger qubit registers. Cannot combine qubit (registers) with gates
        """
        assert self.type==other.type, "Cannot tensor a Gate with a Qubit register"
        if (self.type=="Gate") and (other.type=="Gate"):
            return Gate(tensor(self.array,other.array))
        elif (self.type=="Qubit") and (other.type=="Qubit"):
            return Qubit(tensor(self.array,other.array))

    def ret(self):
        """Returns array of contents of quantum object. Not 'realistic', as
        quantum states cannot be fully observed without collapsing qubit register
        to 1 state. For that, use Qubit.measure()
        """
        return self.array

    def ret_mod(self):
        """Returns array of probabilities of quantum states. Not 'realistic', as
        quantum states cannot be fully observed without collapsing qubit register
        to 1 state. For that, use Qubit.measure()
        """
        return np.abs(np.square(self.array))


    def Generalise(self,n):
        s = self.array
        for i in range(n-1):
            s = tensor(self.array,s)
        self.array = s
#------------------------------------------------------
#---some gates and qubits under QMatrix parent class.
#------------------------------------------------------

# Single Qubit Gates


 
    
        
class Hadamard(QMatrix):
    """Hadamard gate. Takes 0 or 1 state qubit and sets it to equal probability
    superposition.
    """
    def __init__(self,n=1):
        """Initialisation method. n defines number of qubits to act on. Alternatively single qubit
        Hadamards can be tensored together
        """
        QMatrix.__init__(self,"Gate")
        h = np.array([[1,1],[1,-1]])
        hn = h
        for i in range(n-1):
            hn = tensor(h,hn)
        hn = hn*(2**(-0.5*n))
        self.array = np.array(hn)

class V(QMatrix):
    """V gate. Special case of the Phase gate, with phase=pi/2
    """
    def __init__(self):
        """Initialisation method
        """
        QMatrix.__init__(self,"Gate")
        self.array = np.array([[1,0],[0,1j]])

class Phase(QMatrix):
    """Phase gate. Applies complex phase shift to qubit. 
    """
    def __init__(self,phase,n=1):
        """Initialisation method. Phase input defines 
        size of phase shift, n defines number of qubits to act on.
        """
        QMatrix.__init__(self,"Gate")
        
        ph = np.array([[1,0],[0,np.exp(1j*phase)]])
        phn = ph
        for i in range(n-1):
            ph = tensor(ph,phn)
        self.array = ph

class Identity(QMatrix):
    """Identity gate. Leaves qubit register unchanged. Use to represent 'empty' wires
    in quantum circuit diagrams. Typically used by tensoring to other gates
    """
    def __init__(self,n=1):
        """Initialisation method. n defines number of qubits to act on. Alternatively single qubit
        Identity gates can be tensored together
        """
        QMatrix.__init__(self,"Gate")
        self.array = np.identity(2**n)

class PauliZ(QMatrix):
    """Pauli Z gate. Rotates qubit register pi radians about the Z axis of the bloch sphere.
    Special case of Phase shift gate, with phase=pi
    """
    def __init__(self):
        """Initialisation method
        """
        QMatrix.__init__(self,"Gate")
        self.array = np.array([[1,0],[0,-1]])

class PauliY(QMatrix):
    """Pauli Y gate. Rotates qubit register pi radians about the Y axis of the bloch sphere
    """
    def __init__(self):
        """Initialisation method
        """
        SMatrix.__init__(self,"Gate")
        self.array = np.array([[0,-1j],
                               [1j,0]])


class PauliX(QMatrix):
    """Pauli X gate. Rotates qubit register pi radians about the X axis of the bloch sphere
    Classicaly analogous to NOT gate.
    """
    def __init__(self,n=1):
        """Initialisation method
        """
        QMatrix.__init__(self,"Gate")
        self.array = np.flipud(np.identity(2**n))


##############################################################################################
#Control Gates


class Controlled(QMatrix):
    """Generalised controlled gate. Generates controlled versions of any 1 qubit gate,
    where other_gate is another quantum object with type="Gate". 
    """
    def __init__(self,other_gate,n=2):
        """Initialisation method. n defines the number of
        total qubits for the gate (n-1 control qubits + 1 target qubit)
        """
        QMatrix.__init__(self,"Gate")
        self.array = np.identity(2**n)
        self.array[2**n-2,2**n-2] = other_gate.array[0,0]
        self.array[2**n-1,2**n-1] = other_gate.array[1,1]
        self.array[2**n-1,2**n-2] = other_gate.array[1,0]
        self.array[2**n-2,2**n-1] = other_gate.array[0,1]
        


class CNot(QMatrix):
    """Controlled not gate. Equivalent to Controlled(PauliX()).
    """
    def __init__(self,n=2):
        """Initialisation method. n defines the number of
        total qubits for the gate (n-1 control qubits + 1 target qubit)
        """
        QMatrix.__init__(self,"Gate")
        self.array = np.identity(2**n)
        self.array[2**n-2,2**n-2] = 0
        self.array[2**n-1,2**n-1] = 0
        self.array[2**n-1,2**n-2] = 1
        self.array[2**n-2,2**n-1] = 1
        
class Toffoli(QMatrix):
    """Controlled CNot gate. Equivalent to Controlled(PauliX(),3)
    """
    def __init__(self):
        """Initialisation methods
        """
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
    """Controlled phase shift gate. Equivalent to Controlled(Phase()).
    Used a lot for the Quantum Fourier Transform in Shors algorithm
    """
    def __init__(self,phase,n=2):
        """Initialisation method. phase defines phase shift on phase gate, 
        defines the number of total qubits for the gate (n-1 control qubits 
        + 1 target qubit)
        """
        QMatrix.__init__(self,"Gate")
        self.array = np.identity(2**n,dtype=complex)
        self.array[2**n-1,2**n-1]=np.exp(1j*phase)


################################################################
#Other useful gates



class Swap(QMatrix):
    """Swaps the contents of any 2 qubits. Works on entangled qubit registers
    without collapsing them. Mainly useful for shifting qubits around to perform
    controlled gates
    """
    def __init__(self,n=2,index1=0,index2=1):
        """Initialisation method. n defines the total number of qubits being acted on,
        index1 and index2 define the 2 qubits to be swapped
        """
        QMatrix.__init__(self,"Gate")
        self.array = perm_matrix(n,index1,index2)




class Diffusion(QMatrix):
    """Hard coded diffusion gates for Grovers
    """
    def __init__(self,n):
        QMatrix.__init__(self,"Gate")
        N = 2**n
        c = 2.0/N
        self.array = np.full((N, N), c) - np.identity(N)



class Oracle(QMatrix):
    """Hard coded Oracle gate for Grovers
    """
    def __init__(self,reg_size,target):
        QMatrix.__init__(self,"Gate")
        self.array = np.identity(2**reg_size)
        self.array[target,target] = -1

class Noisy(QMatrix):
    """Returns a copy of a Gate, with random noise added. An attempt to more accurately
    simulate real quantum computers, as noise issues are a big deal. However this doesn't
    accurately model how the noise on a real quantum computer emerges, this is just a rough
    attempt. With high amounts of noise, this can break the hermitian and unitary properties
    of gates, and can lead to unpredictable results
    """
    def __init__(self, matrix,level):
        """Initialisation method. other_gate is any other quantum gate. a is a float
        between 0 and 1. 0 -> no noise, 1 -> all noise
        """
        QMatrix.__init__(self,"Gate")
        realnoise = np.random.rand(matrix.array.shape[0],matrix.array.shape[1])
        complexnoise = np.random.rand(matrix.array.shape[0],matrix.array.shape[1])*1j
        noise = realnoise + complexnoise
        self.array = matrix.array*(1-level) + noise*level





class Gate(QMatrix):
    """Generic gate class. Use for defining your own gates. Used by functions that return gates,
    such as multiplication or tensor products.
    """
    def __init__(self,data):
        """Initialise generic gate. data is a matrix. Although no checks are performed,
        matrices must be: Unitary, Hermitian, and size 2^n * 2^n.
        If an input matrix does not satisfy these constraints, unpredictable and wrong 
        things may happen.
        """
        QMatrix.__init__(self,"Gate")
        assert (len(data[0])&(len(data[0])-1)==0) and (len(data[1])&(len(data[1])-1)==0) and (len(data[0])==len(data[1])),"Gate must be square matrix of size 2**n"
        self.array = np.array(data)
      
class Qubit(QMatrix):
    """Generic Qubit class. Used to define qubits. Used as outputs to function such as gate*qubit
    multiplication.
    """
    def __init__(self,data,fock=0):
        """Initialise a qubit or qubit register. If data is an array, set qubit register to that.
        if data is an int, use it as size of qubit register and intialise with optional fock space input
        """
        QMatrix.__init__(self,"Qubit")
        if type(data) is int:
            #If only 1 number is input for data, treat as length of qubit register
            self.array = np.zeros(2**data)
            self.array[fock] = 1
        else:    
            assert (len(data)&(len(data)-1)==0),"Qubit register length must be a power of 2"
            self.array = np.array(data)
        

    def normalise(self):
        """Renormalise qubit register such that probabilities sum to 1
        """
        div = np.sqrt(np.sum(np.square(self.array)))
        a = np.empty(len(self.array))
        a.fill(div)
        self.array = np.divide(self.array,a)

    def measure(self):
        """Measure qubit register. Collapses to 1 definite state. Simulates real measurement, as 
        intermediate values of qubit registers during computation remain unknown.
        """
        pos = np.arange(len(self.array))
        probs = np.abs(np.square(self.array))
        #If probs is not normalised (usually due to rounding errors), re-normalise
        probs = probs/np.sum(probs)
        dist = stats.rv_discrete(values=(pos,probs))
        self.array = np.zeros(self.array.shape)
        self.array[dist.rvs()] = 1
        return self.array

    def measure_cheat(self):
        #Method to measure but ignore state 0. Used to help debug Shors
        pos = np.arange(len(self.array))
        probs = np.abs(np.square(self.array))
        probs[0] = 0
        #If probs is not normalised (usually due to rounding errors), re-normalise
        probs = probs/np.sum(probs)
        dist = stats.rv_discrete(values=(pos,probs))
        self.array = np.zeros(self.array.shape)
        self.array[dist.rvs()] = 1
        return self.array

    def split_register(self):
        """Splits entangle qubit register into individual qubit states. Only works after measurement
        has been acted on a qubit register. Returns a binary string representing the states of each qubit
        """

        outs = np.arange(0,len(self.array),1)
        res = np.array(np.sum(outs*self.array.astype(int)))
        return np.binary_repr(res)

