"""
Input/Output methods
"""

import numpy as np
import math as m
import matplotlib.pyplot as pl

def Hist(qreg):
    x = range(len(qreg))
    pl.bar(x,qreg)
    pl.show()
