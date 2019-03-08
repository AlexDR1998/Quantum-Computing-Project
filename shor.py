import numpy as np
#from gatec import *
from sparse import *
import InOut as IO
import time
import matplotlib.pyplot as plt
from fractions import Fraction
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
	return CPhase(-np.pi*2/(2.0**n))


def H(n=1):
	return Hadamard(n)

def I(n=1):
	return Identity(n)




def HardQFT3():
	#Hard coded 3 qubit qft for testing against
	return (H()&I(2))*(R(2)&I())*(I()&H()&I())*(Swap(3,0,1)*(I()&R(3))*Swap(3,0,1))*(I()&R(2))*(I(2)&H())*Swap(3,0,2)





def HardQFT2():
	#Hard coded 2 Qubit qft
	return (H()&I())*R(2)*(I()&H())*Swap(2,0,1)


def QFT(N):



	def _QFT(n):
		#Method to return gate for n qubit fourier transform
		#Probably best to do recursively
		#print(n)
		if n==2:
			#return ((Identity()&Hadamard())*(R(2))*(Hadamard()&Identity()))
			return ((H()&I())*(R(2))*(I()&H()))#*Swap(2,0,1))

		else:
			g1 = _QFT(n-1)&I()
			g2 = I(n)
			#print(range(n-2))
			
			for i in range(n-2):
				#print(n-i)
				g2 = g2*Swap(n,i,n-2)*(I(n-2)&R(n-i))*Swap(n,i,n-2)
			#IO.Display(g2)
			g3 = (I(n-2)&R(2))*(I(n-1)&H())
			return (g1*g2*g3)
	return _QFT(N)#*Flip(N)


def iQFT(n):
	def _iQFT(n):
		#Inverse fourier transform
		if n==2:
			#return ((Hadamard()&Identity())*(R(2))*(Identity()&Hadamard()))
			return ((I()&H())*(R(2))*(H()&I()))
		else:
			g1 = _iQFT(n-1)&I()
			g2 = I(n)
			for i in range(n-2):
				#print(n-i)
				g2 = Swap(n,i,n-2)*(I(n-2)&R(n-i))*Swap(n,i,n-2)*g2

			g3 = (I(n-1)&H())*(I(n-2)&R(2))
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


def modexp(b,N,qubits):
	#Returns normalised array with modular exponential applied across it.
	def _modexp(b,e,N):
		t = 1
		for i in range(e):
			t = b*t
		return t%N
	
	xs = np.arange(1,2**qubits+1,1)
	g = np.vectorize(lambda x:_modexp(b,x,N))
	f = lambda x:(g(x)-np.mean(x))
	#plt.plot(xs,f(xs))
	#plt.show()
	return f(xs)/np.sqrt(np.sum(np.square(xs).astype(float)))

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





def shor(N,qubits):
	assert N%2!=0, "N must be odd"

	guess = np.random.randint(2,(N-1))
	p = 1
	counter = 0
	while p%2==1 or (guess**(p/2))%N != (N-1):
		counter = counter +1
		guess = np.random.randint(2,(N-1))
		
		#print(guess)
		divisor = GCD(N,guess)
		#print("Does "+str(guess)+" factor "+str(N))

		if divisor!=1:
			#If gcd(N,guess) is not 1, then guess is a factor of N. Lucky
			print("Divisor "+str(divisor)+" factorises "+str(N)+", not found be QFT")
			
			#return guess

		#Run modular exponentiation with base = guess, mod = N
		#p = 0
		#while p==0:
		qreg = Qubit(modexp(guess,N,qubits))
		#print(np.sum(qreg.ret_mod()))
		#Apply QFT to qubit register
		qreg = QFT(qubits)*qreg
		#Measure qubit register to estimate frequency
		#IO.Hist(qreg)
		qreg.measure()

		freq = int(str(qreg.split_register()),2)
		
		#Transform frequency to nearest discrete value
		sample_freq = Fraction(freq,2**qubits)
		p = sample_freq.limit_denominator(N)
		#print(p)
		p = p.denominator
		#print(p)
		
		print("Period guess "+str(p))
		if p%2==0 and ((guess**p/2)%N!=(N-1) and (guess**p/2)%N!=1):
			g1 = GCD(guess**p/2+1,N)
			g2 = GCD(guess**p/2-1,N)
			print(g1)
			print(g2)

			if g1!=1:
				print(str(g1)+" factorises "+str(N))
				#return
			if g2!=1:
				print(str(g2)+" factorises "+str(N))
			if g1!=1 or g2!=1:
				print(str(counter)+" steps")
				return counter
				#return [g1,g2]

	

	#print("Period = "+str(p))

	#print(p)
	
	#if (guess**(p))%N == 1:
	#IO.Hist(qreg)

	




	#print(guess)
	#print(GCD(N,guess))






	#print(extendedGCD(N,guess))

	#print(continued_fraction(guess,N,10))

def main():



	shor(481,9)
	#a = Fraction(1.234567).limit_denominator(1000)

	#print(a)

	#print(extendedGCD(416,93))
	#print(continued_fraction(1,93,416))
	#q1 = Qubit(8)
	#print(q1.ret_mod())
	#q = (Hadamard(3)&Identity(5))*q1
	
	"""
	ft = QFT(10)
	
	qreg = modexp(2311,9123,10)
	q = Qubit(qreg)
	IO.Display(ft)
	IO.Hist(ft*q)
	IO.Hist(ft*Flip(10)*q)
	"""

main()
