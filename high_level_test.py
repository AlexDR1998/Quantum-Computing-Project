import numpy as np
#from gatec import *
from sparse import *
import InOut as IO
import time


import matplotlib.pyplot as plt






def main():
	t1 = time.time()
	#h = Hadamard(1)
	#g1 = Identity()&Hadamard(3)&Identity(8)
	#g2 = Phase(1.2,6)&Hadamard(4)&Identity(2)
	#g3 = Hadamard(2)&Identity(4)&Phase(-0.5,3)&Toffoli()

	#g = g1*g2*g3
	#q1 = Qubit([0,1])
	#q = Hadamard()#Qubit([1,0])
	#q = q1&(Hadamard()*q1)&q1
	#print(q.measure())
	#print(q.split_register())




	o = Oracle(3,1)
	print(o)





	#h = Hadamard(15)
	#p = Phase(1.2,10)
	#g = p*i
	#g = (Identity(10)&Hadamard(1))*p
	#g = (Identity(7)&Hadamard(1)&Identity(7))*p
	#q = Qubit([1,0])
	#q10 = q&q&q&q&q&q&q&q&q&q&q&q&q
	#q2 = h*q10
	t2 = time.time()
	print(t2-t1)
	#p1 = Phase(1,1)
	#.Display(g)
	#IO.Display(g)
	#c2 = CPhase(np.pi/2)
	#c1 = CPhase(np.pi)
	#q = Qubit([1,0])
	#i = Identity()
	#v = V()
	#s = Swap()
	#t = Toffoli()
	#q1 = h*(q&q&q&q)

	#q1.measure()
	#print(q1)
	
	#qft1 = ((i&i&h)*(i&s)*(c2&i)*(i&s)*(c1&i)*(i&h&i)*(c1&i)*(h&i&i))
	#qft = qft1&qft1
	#IO.Display(h&qft1)
	#q_in = q&q&q&q&q&q&q&q&q
	#print(t*q_in)
	#q_in.measure()
	#print(i&h)
	#hist = np.zeros(512)
	#n = 1000
	#for x in range(n):
	#	q_out = qft*q_in
	#print(q_out)
	#	hist = hist + q_out.measure() 
		#print(q_out.measure())

	
	#plt.bar(np.arange(len(hist)),hist/n)
	#plt.show()
	#print(q_out)


	#print(v*h*q)


main()