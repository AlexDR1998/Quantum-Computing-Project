import numpy as np
import scipy as scp
import time
import cmath

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
	



#def main():


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