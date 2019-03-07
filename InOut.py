"""
Input/Output methods
"""

import numpy as np
import math as m
import matplotlib.pyplot as plt
import matplotlib
import random as ran
from gatec import *

#================================   INPUT   ===================================#

def start():
    print('The following methods available for Grovers are:')
    print('input - Enter the number of qubits and target value')
    print('random - Enter the number of qubits and generate a random fock value')
    print('test12 - auto start with 12 qubits for fock 12')
    print('test11 - auto start with 11 qubits for fock 11')
    print('test4 - auto start with 4 qubits for fock 4')
    method = input('Which method would you like: ')
    if method == 'input':
        io = enterVal()
    elif method == 'random':
        io = randVal()
    elif method == 'test12':
        io = [12, 12]
    elif method == 'test11':
        io = [11, 11]
    elif method == 'test4':
        io = [4, 4]
    else:
        print('Not a valid option')

    return io

def enterVal():
    n = int(input('\nHow many qubits? '))
    assert type(n) == int, "n must be an integer greater than or equal to 2"
    assert n >= 2, "n must be an integer greater than or equal to 2"
    N = 2**n
    target = int(input('What Fock space value would you like to find? '))
    assert type(target) == int, "Target must be an integer greater than or equal to 2"
    assert target >= 0, "Target must be an integer greater than or equal to 0"
    assert target <= N-1, "Target must be an integer less than or equal to " + str(N-1)

    io = [n, target]

    return io

def randVal():
    n = int(input('\nHow many qubits? '))
    assert type(n) == int, "n must be an integer greater than or equal to 2"
    assert n >= 2, "n must be an integer greater than or equal to 2"
    N = 2**n
    target = ran.randint(0, N-1)
    io = [n, target]

    return io

#===============================   OUTPUT   ===================================#

def Hist(qreg):
    x = range(len(qreg.ret_mod()))
    plt.bar(x,qreg.ret_mod())
    plt.show()


def Graph(qreg):

	x = range(len(qreg.ret_mod()))
	plt.plot(x,qreg.ret_mod())
	plt.show()


def Display(Gate):
	#Function to plot gate matrix as image.
	m = Gate.ret()
	plt.imshow(complex_array_to_rgb(m),cmap="nipy_spectral")
	plt.show()


def complex_array_to_rgb(X, theme='dark', rmax=None):
	#maps array of complex numbers to colours. Taken from stack overflow:
	#https://stackoverflow.com/questions/15207255/is-there-any-way-to-use-bivariate-colormaps-in-matplotlib

	absmax = rmax or np.abs(X).max()
	Y = np.zeros(X.shape + (3,), dtype='float')
	Y[..., 0] = np.angle(X)/(2 * np.pi) % 1
	if theme == 'light':
		Y[..., 1] = np.clip(np.abs(X) / absmax, 0, 1)
		Y[..., 2] = 1
	elif theme == 'dark':
		Y[..., 1] = 1
		Y[..., 2] = np.clip(np.abs(X) / absmax, 0, 1)
	Y = matplotlib.colors.hsv_to_rgb(Y)
	return Y
