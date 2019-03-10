import numpy as np
from gatec import *
#from sparse import *
import InOut as IO
import time
import matplotlib.pyplot as plt
from fractions import Fraction


#Code to implement Shor's algorithm for polynomial time factoring of numbers.
#Broken up into various subroutines (some quantum some classical)

#----------------------------------------------------
#--- Quantum subroutines --------------- AR ---------
#----------------------------------------------------

#--- Helper functions for QFT ---

def Flip(n,noise=0):
	#Applies lots of Swap gates to flip qubit register
	if noise==0:
		g = Identity(n)
		for i in range(n//2):
			if i!=(n-1-i):
				g = g*Swap(n,i,n-1-i)
		return g
	else:
		return Noisy(Flip(n,0),noise)

#Simple gates, reduced to 1 letter function calls for clarity
def R(n,noise=0.1):
	if noise==0:
		return CPhase(-np.pi*2/(2.0**n))
	else:
		return Noisy(CPhase(-np.pi*2/(2.0**n)),noise)

def H(n=1,noise=0.1):
	if noise==0:
		return Hadamard(n)
	else:
		return Noisy(Hadamard(n),noise)
def I(n=1,noise=0):
	if noise==0:
		return Identity(n)
	else:
		return Noisy(Identity(n),noise)

#--- QFT defined recursively ---

"""
def QFT(N):

	def _QFT(n):
		#Method to return gate for n qubit fourier transform
		if n==2:
			#Base case
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
	#Some ambiguity whether Flip gate is needed here. Appears to work with or without
	return _QFT(N)*Flip(N)

"""

def QFT(N):

	def _QFT(n):
		#Method to return gate for n qubit fourier transform
		if n==2:
			#Base case
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
	#Some ambiguity whether Flip gate is needed here. Appears to work with or without
	return _QFT(N)*Flip(N)



#--- Inverse QFT defined recursively ---

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


"""

def HardQFT3():
	#Hard coded 3 qubit qft for testing against
	return (H()&I(2))*(R(2)&I())*(I()&H()&I())*(Swap(3,0,1)*(I()&R(3))*Swap(3,0,1))*(I()&R(2))*(I(2)&H())*Swap(3,0,2)


def HardQFT2():
	#Hard coded 2 Qubit qft
	return (H()&I())*R(2)*(I()&H())*Swap(2,0,1)

"""




#-------------------------------------------------------------------------
#--- Classical subroutines -----------------------------------------------
#-------------------------------------------------------------------------



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

def isPrime(n):
	#Method to check if a number is prime. Used to check input to shor's is not prime
	if n <= 3:
		return n>1
	elif (n%2==0) or (n%3==0):
		return False
	else:
		i = 5
		while i*i <= n:
			if (n%i==0) or (n%(i+2)==0):
				return False
			i = i+6
		return True


def modexp(b,N,qubits):
	#Returns normalised array with modular exponential applied across it.
	def _modexp(b,e,N):
		t = 1
		for i in range(e):
			t = b*t
		return t%N
	
	xs = np.arange(0,2**qubits,1)
	g = np.vectorize(lambda x:_modexp(b,x,N))
	f = lambda x:(g(x))#-np.mean(x))
	#plt.plot(xs,f(xs))
	#plt.show()
	#print (np.sum(np.square(np.abs(f(xs)/np.sqrt(np.sum(np.square(xs).astype(float)))))))
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
	assert isPrime(N)==False, "Cannot factorise prime number"+str(N)

	guess = np.random.randint(2,(N-1))
	p = 1
	counter = 0
	found = False

	qubit_profile = np.zeros(2**qubits)

	
	while found==False:#p%2==1 or (guess**(p/2))%N != (N-1):
		counter = counter +1
		guess = np.random.randint(2,(N-1))
		
		#print(guess)
		divisor = GCD(N,guess)
		#print("Does "+str(guess)+" factor "+str(N))

		if divisor!=1:
			
			#If gcd(N,guess) is not 1, then guess is a factor of N. Lucky
			#print("Divisor "+str(divisor)+" factorises "+str(N)+", not found by QFT")
			#found = True
			pass
			#return [divisor,divisor,counter]

		#Run modular exponentiation with base = guess, mod = N
		qreg = Qubit(modexp(guess,N,qubits))
		#print(np.sum(qreg.ret_mod()))
		#Apply QFT to qubit register
		qreg = iQFT(qubits)*qreg
		
		#qubit_profile = np.vstack((qubit_profile,qreg.ret_mod()))
		#Measure qubit register to estimate frequency

		#IO.Hist(qreg)
		qreg.measure_cheat()
		freq = int(str(qreg.split_register()),2)
		
		#Transform frequency to nearest discrete value
		sample_freq = Fraction(freq,2**qubits)
		p = sample_freq.limit_denominator(N).denominator
		#print("Period guess "+str(p))

		#if period is odd or is trivial square root, start again

		if p%2==0 and ((guess**p/2)%N!=(N-1) and (guess**p/2)%N!=1):
			g1 = GCD(guess**p/2+1,N)
			g2 = GCD(guess**p/2-1,N)
			#print(g1)
			#print(g2)

			if g1!=1:
				#print(str(g1)+" factorises "+str(N))
				found = True
				#return
			if g2!=1:
				#print(str(g2)+" factorises "+str(N))
				found = True
			if g1!=1 or g2!=1:
				#print(str(counter)+" steps")
				#return counter
				#plt.matshow(qubit_profile)
				#plt.show()
				return [g1,g2,counter]
	return shor(N,qubits)
	

	



#---------------------------------------------------------
#--- Methods for testing shors ---------------------------
#---------------------------------------------------------

def coPrimes(bits,bits_lower = 1):
	#Randomly generates a product of 2 primes, n, such that
	# n <2**bits
	#Also returns the 2 prime factors

	#print(list(filter(isPrime,range(2**bits))))
	primes = np.array(list(filter(isPrime,range(2**bits))))[1:]
	prime_products = np.outer(primes,primes)

	res = 0
	while (res >2**bits) or (res < 2**(bits_lower)):
		x = np.random.randint(0,len(primes))
		if x==0:
			y=0
		else:
			y = np.random.randint(0,x+1)
		res = prime_products[x,y]
	return res,primes[x],primes[y]


def step_test(lower,upper,its=10):
	#Method for testing shors on randomly generated prime products.
	#Plots number of QFTs applied against size of Qubit register
	hits = 0
	misses = 0
	
	steps = np.zeros((its,len(range(lower,upper))))

	for bits in (range(lower,upper)):
		print("Runnint tests for "+str(bits)+" Qubits")
		for x in range(its):
			print(x)
			coprime,prime1,prime2 = coPrimes(bits+1,bits-2)
			res = shor(coprime,bits)
			steps[x,bits-lower] = res[2]
			if (res[0]==prime1) or (res[0]==prime2) or (res[1]==prime1) or (res[1]==prime2):
				hits = hits+1
			else:
				misses = misses +1
	ms = np.mean(steps,axis=0)
	ers = np.std(steps,axis=0)
	plt.errorbar(x=range(lower,upper),y=ms,yerr=ers,fmt='o')
	plt.show()
	print(str(hits)+" succesful factorings")
	print(str(misses)+" failures")



def main():

	#step_test(4,11,20)
	#print(coPrimes(7,4))
	#print(shor(7*17,8))
	#print(isPrime(209))
	#for x in range(10):
	#	print(coPrimes(8))
	#ls = filter(isPrime,range(1000))
	#print(ls)
	#a = Fraction(1.234567).limit_denominator(1000)

	#print(a)

	#print(extendedGCD(416,93))
	#print(continued_fraction(1,93,416))
	#q1 = Qubit(8)
	#print(q1.ret_mod())
	#q = (Hadamard(3)&Identity(5))*q1
	
	
	#ft = QFT(6)
	
	#qreg = modexp(2311,9123,10)
	#q = Qubit(qreg)
	IO.Display(QFT(7))#*QFT(7)*QFT(7)*QFT(7))
	#IO.Hist(ft*q)
	#IO.Hist(ft*Flip(10)*q)
	

main()
