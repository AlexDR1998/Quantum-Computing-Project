import numpy as np

from sparse import *
import InOut as IO
import time











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
	print(extendedGCD(416,93))
	print(continued_fraction(1,93,416))

main()
