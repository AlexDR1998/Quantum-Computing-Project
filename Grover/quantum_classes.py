import numpy as np

from binar_algebra_fct import (
	bin2int, int2bin
)


class Qubit:
	"""
	Class of a pure quantum bit
	It includes two representation : 
	The ket representation : |01001>
	The vector space : 0 x 1 x 0 x 0 x 1 (a vector of length 2^5 with only one 1)
	"""
	
	def __init__(self,type='k',var=[0]): 
		"""type is 'k' for ket and 'v'for vector. ket should be given as an array of size n, and vector an array of size 2**n"""
	
		if type=='k':
			self.ket=array(var)
			self.vec=zeros(2**(len(self.ket)),int) # the size of the space is 2**n , where n=len(k)
			self._vec_coord=bin2int(str(array2string(self.ket, precision=0, separator='')[1:-1])) # We find which coord has to be 1
			self.vec[self._vec_coord]+=1 # we update it
		if type=='v':
			self.vec=array(var)
			self._lenket=int(log2(len(self.vec))) # the size of the ket is log2(vec)
			self.ket=array(list(int2bin(argmax(self.vec),self._lenket)),int)
			
			
class EQubits:
	"""
	Class corresponding to entangled states of Qubits
	"""
	
	def __init__(self,type='k',var=[[0,1]]):
		"""
		type is 'k' for ket and 'v'for vector. 
		ket should be given as [[qubit1.ket,p1],[qubit2.ket,p2]] and sum(pi^2)=1, 
		vector a normalized array of size 2**n
		"""
		
		if type=='k':
			#self.ket=array([[ ]])
			self.ket=empty((1,2))
			self.vec = zeros(2**(len(var[0][0])))
			for i in var:
				self._tempqb=(Qubit('k',(i[0])))
				self.ket = append(self.ket,array([( self._tempqb.ket,i[1])],object),axis=0)
				self.vec+=(self._tempqb.vec)*i[1]
			self.ket=self.ket[1:]
			#self.ket=array([(Qbit1,p1),...,(QbitN,pN)])
			#self.vec=array([0,0,pk,0,0,...,pj,0,...,0])
			
		if type=='v':
			self.vec=array(var)
			self._lenket=int(log2(len(self.vec))) # the size of the ket is log2(vec)
			self.ket=empty((1,2))
			for i in flatnonzero(self.vec):
				self._zero=zeros(len(self.vec))
				self._zero[i]=1
				self._qtemp=Qubit('v',self._zero)
				
				#print array([[self._qtemp.ket,self.vec[i]]])
				self.ket = append(self.ket,array([(self._qtemp.ket,self.vec[i])],object),axis=0)
			self.ket=self.ket[1:]
			#self.ket = array([Qbit1,p1],...,[QbitN,pN])
	
	def update_vec(self):
		"Update the QUbit when you send a new vector) - TODO"

	def number_states(self):
		"Returns the number of entangled vector in the Qubit"
		return len(self.ket)
		
def hadamard(n):
	"""Returns the n-th Hadamard matrix """
	
	H1=1/np.sqrt(2)*np.array([1,1,1,-1]).reshape((2,2))
	Hn=H1
	for i in range(n-1):
		Hn=np.kron(Hn,H1)
	return Hn