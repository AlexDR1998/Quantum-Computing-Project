import numpy as np
import scipy as scp
import time
import cmath
from lazyarray import larray
#from lazyarray import matmul
"""
A file for storing various mathematical helper functions that could be used in various places.
 - AR
"""
def lazy_mul_gate(b,a):

	#Dimension of output
	a0 = a.shape[0]
	a1 = a.shape[1]
	b0 = b.shape[0]
	b1 = b.shape[1]
	outdim = (b0,a1)

	
	#Calculate output matrix
	def mul(i,j):
		n=0
		element = 0
		for n in range( outdim[0]):
			element += (a[i,n]*b[n,j])
		return element	
	
	output = larray(mul,shape=outdim)
	
	if output.shape[0]==1:
		output = output[0]
	return(output)

def lazy_mul(b,a):
	#Dimension of output
	a0 = a.shape[0]
	a1 = a.shape[1]
	b0 = b.shape[0]
	b1 = b.shape[1]
	outdim = (b0,a1)
	print(b.evaluate())
	print(a.evaluate())
	print(outdim)
	#output = np.zeros(outdim)
	#print("gate"+str(b.evaluate()))
	#print("qubit"+str(a.evaluate()))
	#print(outdim)

	def mul(i,j):
		def k(i,j):
			return j
		lis = larray(k,shape = (1,b1))

		elem = sum(map(lambda n: b[i][n]*a[n],lis[0]))
		#print(elem)
		return elem
	
		####
		#elem = 0
		#for n in range(a.shape[0]):
		#	elem += b[j][n]*a[n]
		#print(elem)
		#return elem
		#####

	output = larray(mul,shape = outdim)
	#print("HEY")
	#output = larray([output])
	#print (output.evaluate())
	return output


def tensor(b,a):

	"""
	Function that takes the tensor product of 2 matrices or vectors.
	Behaves identically to np.kron()
	"""
	#Dimension of output
	
	a0 = a.shape[0]
	a1 = a.shape[1]
	b0 = b.shape[0]
	b1 = b.shape[1]
	outdim = (a0*b0,a1*b1)

	#Initialise output matrix with zeros
	output = np.zeros(outdim,dtype=complex)

	#Calculate output matrix
	for x in range(outdim[0]):
		for y in range(outdim[1]):
			output[x,y] = a[x%a0,y%a1]*b[x//a0,y//a1]

	return output
	#If output matrix is (1,n), then just convert to n vector

def tensor_lazy(b,a):

	"""
	Function that takes the tensor product of 2 matrices or vectors.
	Behaves identically to np.kron()
	"""
	#Convert any input vectors to (1,n) matrices
	#Dimension of output

	a0 = a.shape[0]
	a1 = a.shape[1]
	b0 = b.shape[0]
	b1 = b.shape[1]
	outdim = (a0*b0,a1*b1)


	#Calculate output matrix
	def kron(i,j):
		return a[i%a0,j%a1]*b[i//a0,j//a1]

	output = larray(kron,shape=outdim)
	return output

def tensor_lazied(b,a):

	"""
	Function that takes the tensor product of 2 matrices or vectors.
	Behaves identically to np.kron()
	"""
	#Convert any input vectors to (1,n) matrices
	#Dimension of output

	a0 = a.shape[0]
	a1 = a.shape[1]
	b0 = b.shape[0]
	b1 = b.shape[1]
	outdim = (a0*b0,a1*b1)


	#Calculate output matrix
	def kron(i,j):
		return a[i%a0,j%a1]*b[i//a0,j//a1]

	output = larray(kron,shape=outdim)
	return output.evaluate()


def tensor_sparse_gate(A,B):
	from scipy import sparse as sp
	#return sp.kron(A,B)
	#sp.kron and tensor_sparse give the same result

	# B is fairly dense, use BSR
	A = sp.csc_matrix(A,copy=True)

	output_shape = (A.shape[0]*B.shape[0], A.shape[1]*B.shape[1])

	if A.nnz == 0 or B.nnz == 0:
	    # kronecker product is the zero matrix
	    return sp.coo_matrix(output_shape)

	B = B.toarray()
	data = A.data.repeat(B.size).reshape(-1,B.shape[0],B.shape[1])
	data = data * B

	return sp.bsr_matrix((data,A.indices,A.indptr), shape=output_shape)


def tensor_sparse_qubit(A,B):

	#For sparse qubits, probably just easier to convert to dense arrays.
	#As 1D, gains made from sparse optimisations are minimal

	a = A.toarray()[0]
	b = B.toarray()[0]

	#Dimension of output
	a0 = a.shape[0]
	b0 = b.shape[0]
	outdim = (1,a0*b0)

	#Initialise output matrix with zeros
	output = np.zeros(outdim,dtype=complex)

	#Calculate output matrix
	for x in range(outdim[0]):
		output[0,x] = a[x%a0]*b[x//a0]

	return(output)




def perm_matrix(n,index1,index2):
	#generates a permutation matrix from a list of pairs of numbers to swap
	

	assert index1!=index2, "Cant swap qubit with itself"
	assert (index1<n) and (index2<n), "Cant swap qubits beyond size of gate"
	index1 = n-index1-1
	index2 = n-index2-1
	b = 2**index1 + 2**index2

	swaps = []
	for x in range(2**n):
		for y in range(x):
			if((x^y==b) and (count_bits(x)==count_bits(y))):
				swaps.append((x,y))

	#print(swaps)
	size = 2**n
	i = np.identity(size)
	for pairs in swaps:
		temp = i[pairs[0]].copy()
		i[pairs[0]] = i[pairs[1]]
		i[pairs[1]] = temp
	return i

def count_bits(n):
	if n==0:
		return 0
	else:
		return (n&1)+count_bits(n>>1)
