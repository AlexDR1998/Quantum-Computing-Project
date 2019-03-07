import numpy as np
from gatec import *
#from sparse import *
import InOut as IO
import time

#Code to implement Shor's algorithm for polynomial time factoring of numbers.
#Broken up into various subroutines (some quantum some classical)




def Flip(n):
	#Applies lots of Swap gates to basically flip qubit register
	g = Identity(n)
	for i in range(n//2):
		if i!=(n-1-i):
			g = g*Swap(n,i,n-1-i)
			#print(i)
	return g

def R(n):
	#Helper function for QFT()
	return CPhase(np.pi*2/(2**n))

def HardQFT3():
	#Hard coded 3 qubit qft for testing against
	I = Identity()
	H = Hadamard()
	return (H&Identity(2))*(R(2)&I)*(I&H&I)*(Swap(3,0,1)*(I&R(3))*Swap(3,0,1))*(I&R(2))*(Identity(2)&H)*Swap(3,0,2)

def HarderQFT3():
	w = np.sqrt(1j)
	w3 = w**3
	w5 = w**5
	w7 = w**7
	return Gate(1/np.sqrt(8)*np.array([[1,1,1,1,1,1,1,1],
				 		   	  [1,w,1j,w3,-1,w5,-1j,w7],
				 		   	  [1,1j,-1,-1j,1,1j,-1,-1j],
						   	  [1,w3,-1j,w,-1,w7,1j,w5],
				 		   	  [1,-1,1,-1,1,-1,1,-1],
				 		   	  [1,w5,1j,w7,-1,w,-1j,w**3],
				 		   	  [1,-1j,-1,1j,1,-1j,-1,1j,1],
				 		   	  [1,w7,-1j,w5,-1,w3,1j,w]]))

def HardQFT2():
	#Hard coded 2 Qubit qft
	I = Identity()
	H = Hadamard()
	return (H&I)*R(2)*(I&H)*Swap(2,0,1)


def QFT(n):

	def _QFT(n):
		#Method to return gate for n qubit fourier transform
		#Probably best to do recursively
		if n==2:
			#return ((Identity()&Hadamard())*(R(2))*(Hadamard()&Identity()))
			return ((Hadamard()&Identity())*(R(2))*(Identity()&Hadamard()))

		else:
			g1 = _QFT(n-1)&Identity()
			g2 = Identity(n)
			for i in range(n-2):
				#print(n-i)
				g2 = g2*Swap(n,i,n-1)*(Identity(n-2)&R(n-i))*Swap(n,i,n-1)
			#IO.Display(g2)
			g3 = (Identity(n-2)&R(2))*(Identity(n-1)&Hadamard())
			return (g1*g2*g3)
	return _QFT(n)*Flip(n)


def iQFT(n):
	def _iQFT(n):
		#Inverse fourier transform
		if n==2:
			#return ((Hadamard()&Identity())*(R(2))*(Identity()&Hadamard()))
			return ((Identity()&Hadamard())*(R(2))*(Hadamard()&Identity()))
		else:
			g1 = _iQFT(n-1)&Identity()
			g2 = Identity(n)
			for i in range(n-2):
				#print(n-i)
				g2 = Swap(n,i,n-1)*(Identity(n-2)&R(n-i))*Swap(n,i,n-1)*g2

			g3 = (Identity(n-1)&Hadamard())*(Identity(n-2)&R(2))
			return (g3*g2*g1)
	return Flip(n)*_iQFT(n)


def GCD(x,y):
	#Returns greatest common divisor of x and y. Pretty fast, shouldn't be the bottleneck anyway
	if x==y:
		return x
	while(y): 
		x,y = y,x%y 
	return x

def extendedGCD(x,y):
	#Returns array representing continued fraction expansion of x/y
	fracs = []
	if x==y:
		return fracs
	while(y):
		fracs.append(x//y)
		x,y = y,x%y
	return np.array(fracs)

#def isPrime(n):
#Method to check if a number is prime. Used to check input to shor's is not prime


def continued_fraction(y, Q, N):
	#Not 100% sure what this does
	fractions = extendedGCD(y, Q)
	depth = 2

	def partial(fractions, depth):
		c = 0
		r = 1
		for i in reversed(range(depth)):
			tR = fractions[i] * r + c
			c = r
			r = tR
		return c
	r = 0
	for d in range(depth, len(fractions) + 1):
		tR = partial(fractions, d)
		if tR == r or tR >= N:
			return r
		r = tR
	return r





def shor(N):
	assert N%2!=0, "N must be odd"
	
	guess = np.random.randint(1,N)
	print(guess)
	divisor = GCD(N,guess)

	if divisor!=1:
		#If gcd(N,guess) is not 1, then guess is a factor of N. Lucky
		return guess


	#print(guess)
	print(GCD(N,guess))






	#print(extendedGCD(N,guess))

	print(continued_fraction(guess,N,10))

def main():
	
	#shor(15)
	#print(extendedGCD(416,93))
	#print(continued_fraction(1,93,416))
	#q1 = Qubit(3)
	#print(q1.ret_mod())
	#q = (Hadamard()&Identity(2))*q1
	#IO.Hist(q)
	
	#print(q1)
	#print(q2)
	ft = QFT(3)
	ift = iQFT(3)
	ft2 = HardQFT3()
	#ift = iQFT(3)
	IO.Display(ft)
	IO.Display(ft2)
	IO.Display(ift)
	IO.Display(ift*ft)
	#IO.Hist(ft*q)
	#IO.Display(Swap(4,0,3)*Swap(4,1,2))#*Swap(5,1,3))
	#IO.Display(Flip(8))
	#IO.Graph(ft*q)

main()
