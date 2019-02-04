"""
Contains all the functions relative to the perceptron
"""

import numpy as np

from binar_algebra_fct import (
	bin2int, int2bin, hamming_distance
)

def perc_spl(a,x,n):	
	"Return 1 if the binary representation of a and x differs by at least n/2"
	
	ab=np.binary_repr(a, width=n)
	xb=np.binary_repr(x, width=n)
	return int(hamming_distance(ab, xb)<=np.log2(n)/2)

def mat_hamming(N):	
	"Generation de la matrice de Hamming Distance for N=2**n bitstrings"
	
	hamming=np.zeros((N,N))	
	for i in range(N):
		for j in range(N):
			ib=np.binary_repr(i, width=N)
			jb=np.binary_repr(j, width=N)
			hamming[i][j]=hamming_distance(ib,jb)
	return hamming
	
def mat_teacher_bin(n,thres=-1):
	"""
	Generation of the membership query matrix. Every column corresponds to the teacher answer 
	of a particular weight distribution (the column), to a given input ( the row)
	"""
	if thres==-1:
		thres=int(np.log2(n))/2
	teachers=np.zeros((n,n),complex)	
	for i in range(n):
		for j in range(n):
			ib=np.binary_repr(i, width=n)
			jb=np.binary_repr(j, width=n)
			teachers[i][j]=(-1)**(hamming_distance(ib,jb)<=thres)
	return teachers	
	
def mat_teacher_phase(n):
	"""
	Generation of the membership query matrix with phase instead of +/-1. Every column corresponds to the teacher answer 
	of a particular weight distribution (the column), to a given input ( the row), depending on the hamming distance
	"""
	
	ni=int(log2(n))
	res=complex(0+1j*0)
	for i in range(n):
		d=hamming[weight][i]
		c=tab_coef[i]
		if tab_teacher[i]:
			res+=complex(exp(1j*(pi*(d)/ni))) #First version : pi/n	
		else :
			res+=complex(exp(1j*(pi*(d/ni+1))))
			
	teachers=np.zeros((n,n))	
	for i in range(n):
		for j in range(n):
			ib=np.binary_repr(i, width=n)
			jb=np.binary_repr(j, width=n)
			teachers[i][j]=(-1)**(hamming_distance(ib,jb)<=int(np.log2(n))/2)
	return teachers	