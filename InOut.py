"""
Input/Output methods
"""

import numpy as np
import math as m
import matplotlib.pyplot as plt
import matplotlib
import random as ran
import sys
from dense import *
#from sparse import *

#================================   INPUT   ===================================#
def vers():
    '''
    Input check for how to run Grover simulation; either testing or standard/noisy Grover
    '''
    print('\nOptions:')
    print('test1 - Gather data and plot runlength versus number of qubits for fixed target value')
    print('test2 - Gather data and graph runlength versus target value for fixed qubit register')
    print('run - start as usual with no noise')
    print('noisy - run with noise')
    check = input(('Enter your choice: '))

    return check

def start():
    '''
    For running Grovers, provides pre set options or allowed to input reg size and fock target
    '''
    print('The following methods available for Grovers are:')
    print('input - Enter the number of qubits and target value(in Fock space)')
    print('random - Enter the number of qubits and generate a random Fock space value')
    print('test12 - auto start with 12 qubits searching for |12>')
    print('test11 - auto start with 11 qubits searching for |11>')
    print('test4 - auto start with 4 qubits searching for |4>')
    method = input('Enter the method would you like: ')
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
        sys.exit()

    return io

def enterVal():
    '''
    Gathers the users chosen reg size and target fock value
    '''
    n = int(input('\nEnter the number of qubits? '))
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
    '''
    Gathers users chosen reg size an creates a random fock target
    '''
    n = int(input('\nEnter the number of qubits? '))
    assert type(n) == int, "n must be an integer greater than or equal to 2"
    assert n >= 2, "n must be an integer greater than or equal to 2"
    N = 2**n
    target = ran.randint(0, N-1)
    io = [n, target]

    return io

def gnoise():
    '''
    Gathers the users chosen noise level
    '''
    print('\nThe smallest noise is 0 and largest is 1')
    noise = float(input('Enter a value for noise: '))
    assert noise >= 0, "noise value must be between 0 and 1"
    assert noise <= 1, "noise value must be between 0 and 1"

    return noise

#===============================   OUTPUT   ===================================#
def printOut(q, target):
    '''
    Prints the output and target of Grovers in both binary and standard
    '''
    print('\nThe state of the ouput(in binary) is |' + str(q.split_register()) + '>')
    print('In Fock space this is |' + str(int(str(q.split_register()), 2)) + '>')
    print('The target state(in binary) was |' + str(bin(target)[2:]) + '>')
    print('In Fock space this is |' + str(target) + '>')

def timeplotn(x, y):
    '''
    Used for plotting Grover test results for varying reg size
    '''
    plt.plot(x, y)
    plt.title('Time for Grovers to run versus the number of qubits')
    plt.xlabel('Number of qubits')
    plt.ylabel('Time, s')
    plt.show()

def timeplottar(x, y):
    '''
    Used for plotting Grover test results for varying fock target
    '''
    plt.plot(x, y)
    plt.title('Time for Grovers to run versus the target Fock value for 10 qubits')
    plt.xlabel('Fock Value')
    plt.ylabel('Time, s')
    plt.show()

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
