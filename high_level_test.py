import numpy as np
#from gatec import *
from sparse import *
import InOut as IO
import time
from fractions import Fraction


import matplotlib.pyplot as plt



def extendedGCD(x,y):
	#Returns array representing continued fraction expansion of x/y
	fracs = []
	if x==y:
		return fracs
	while(y):
		fracs.append(x//y)
		x,y = y,x%y
	return np.array(fracs)



def main():
	
	n = 123
	d = 251
	
	a = Fraction(n,d)
	print(a.limit_denominator(100).numerator)
	
	print(extendedGCD(n,d))

main()