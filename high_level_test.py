import numpy as np
from gatec import *
import matplotlib.pyplot as plt

def main():
	h = Hadamard()
	b = CPhase(np.pi/4)
	q = Qubit([1,0])
	i = Identity()
	v = V()

	#q1 = h*(q&q&q&q)

	#q1.measure()
	#print(q1)
	alg = (h&h&h)*(b&i)*(v&i&v)*(i&h&i)
	
	q_in = q&q&q

	hist = np.zeros(8)
	n = 1000
	for x in range(n):
		q_out = alg*q_in
		hist = hist + q_out.measure() 
		#print(q_out.measure())

	
	plt.bar(np.arange(len(hist)),hist/n)
	plt.show()
	#print(q_out)


	#print(v*h*q)


main()