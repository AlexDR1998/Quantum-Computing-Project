import numpy as np
import scipy as scp
import time
import cmath
from scipy import sparse as sp
from lazyarray import larray

"""
A file for storing various mathematical helper functions that could be used in various places.
 - AR
"""


def tensor(b,a):
	
	"""
	Function that takes the tensor product of 2 matrices or vectors.
	Behaves identically to np.kron()
	"""
	#Convert any input vectors to (1,n) matrices
	if (len(a.shape)==1):
		a = np.array([a])
	if (len(b.shape)==1):
		b = np.array([b])

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
	
	
	#If output matrix is (1,n), then just convert to n vector
	if output.shape[0]==1:
		output = output[0]
	return(output)
	
def tensor_lazy(b,a):
	
	"""
	Function that takes the tensor product of 2 matrices or vectors.
	Behaves identically to np.kron()
	"""
	#Convert any input vectors to (1,n) matrices
	if (len(a.shape)==1):
		a = np.array([a])
	if (len(b.shape)==1):
		b = np.array([b])

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
	
	
	#If output matrix is (1,n), then just convert to n vector
	if output.shape[0]==1:
		output = output[0]
	return(output)
	
def tensor_sparse(A,B):
	#Numpy implementation for now - will replace with my own soon - AR
	#return sp.kron(a,b)

	
	# B is fairly dense, use BSR
	A = sp.csr_matrix(A,copy=True)

	output_shape = (A.shape[0]*B.shape[0], A.shape[1]*B.shape[1])

	if A.nnz == 0 or B.nnz == 0:
	    # kronecker product is the zero matrix
	    return sp.coo_matrix(output_shape)

	B = B.toarray()
	data = A.data.repeat(B.size).reshape(-1,B.shape[0],B.shape[1])
	data = data * B

	return sp.bsr_matrix((data,A.indices,A.indptr), shape=output_shape)

def perm_matrix(n,index1,index2):
	#generates a permutation matrix from a list of pairs of numbers to swap
	assert index1!=index2, "Cant swap qubit with itself"
	assert (index1<n) and (index2<n), "Cant swap qubits beyond size of gate"
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


#def main():

	#q1 = 0
	#q2 = 1

	#b = 2**q1+2**q2

	#b = int("101",2)
	#swaps = []
	#for x in range(8):
	#	for y in range(x):
	#		if((x^y==b) and (count_bits(x)==count_bits(y))):
	#			swaps.append((x,y))
	#			#print(x,y)
	#print(perm_matrix(3,swaps))
	"""
	Stuff to test any functions defined here - DELETE LATER -AR
	"""

	#x = np.array([1,2,3])
	#y = np.array([1,-1])
	#x = np.random.rand(1,5)
	#y = np.random.rand(1,3)
	#t1 = time.time()
	#a = tensor(x,y)
	#t2 = time.time()

	#print("Time taken: "+str(t2-t1))
	#print(x)
	#print(y)
	#print(a)
	#print (a==np.kron(x,y)).all()
	

	#for n in range(10):
	#	x = 2**n
	#	print(x&(x-1))

#main()