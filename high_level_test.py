import numpy as np
from gatec import *

import matplotlib.pyplot as plt

def main():
	h = Hadamard()
	c2 = CPhase(np.pi/2)
	c1 = CPhase(np.pi)
	q = Qubit([1,0])
	i = Identity()
	v = V()
	s = Swap()
	#q1 = h*(q&q&q&q)

	#q1.measure()
	#print(q1)
	
	qft = (i&i&h)*(i&s)*(c2&i)*(i&s)*(c1&i)*(i&h&i)*(c1&i)*(h&i&i)
	

	q_in = q&q&q
	#q_in.measure()

	hist = np.zeros(8)
	n = 1000
	for x in range(n):
		q_out = qft*q_in
	#print(q_out)
		hist = hist + q_out.measure() 
		#print(q_out.measure())

	
	plt.bar(np.arange(len(hist)),hist/n)
	plt.show()
	#print(q_out)


	#print(v*h*q)


main()