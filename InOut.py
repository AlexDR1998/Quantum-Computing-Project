"""
Input/Output methods
"""

import numpy as np
import math as m
import matplotlib.pyplot as plt
import matplotlib
from gatec import *


def Hist(qreg):
    x = range(len(qreg))
    plt.bar(x,qreg)
    plt.show()


def Graph(qreg):
	x = range(len(qreg))
	plt.plot(x,qreg)
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