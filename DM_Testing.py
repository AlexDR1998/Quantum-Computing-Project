import numpy as np
import math as m
from gatec import *

def main():

    pio2 = m.pi/2.0

    q1 = Qubit([1, 0])
    q2 = Qubit([0, 1])
    h = Hadamard()
    I = Identity()
    cNOT = CNot()
    X = PauliX()
    S = Phase(pio2)

    qreg = q1&q1
    h1 = I&h #h on qubit 1
    hh = h&h #h on both quibits
    XX = X&X #PauliX on both
    SS = S&S #Phase on both

    #NOTE: By first qubit I mean the one that has the large solid dot of the
    #      cNOT gate in circuit diagrams, usually seen as the lower "wire"

    #No phase on any = [0,0,0,1] ; this is the 11 state
    #Phase on both should give 00 state
    #Phase on first should equal 01 state
    #Phase on second should equal 10 state

    qf = hh*qreg

    #For 2 quibits in here the phase gate is inserted for either 1, both or neither quibits
    qf = SS*qf #phase both

    #If it is inserted then it has to be mirrored after the cNOT group

    qf = h1*qf
    qf = cNOT*qf
    qf = h1*qf

    #Here is where the mirroring occurs
    qf = SS*qf

    #This section is untouched in the 2 quibit Grover

    qf = hh*qf
    qf = XX*qf

    qf = h1*qf
    qf = cNOT*qf
    qf = h1*qf

    qf = XX*qf
    qf = hh*qf

    qf.measure()
    print(qf)

main()
