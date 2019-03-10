import numpy as np
import math as m
import time as t
import sys
import InOut as IO
import GroverGateWise as G

def main():
    check = IO.vers()

    if check == 'run':
        args = IO.start()
        G.run(args)

    elif check == 'noisy':
        args = IO.start()
        noise = IO.gnoise()
        G.runnoisy(args, noise)

    elif check == 'test':
        # --- Time for n qubits ---
        args = np.zeros(2)
        nlist = [i for i in range(2,14)]
        time1 = np.zeros(len(nlist))
        target = 1
        for i in range(len(nlist)):
            print('Running Grovers for ' + str(nlist[i]) + ' qubits for Fock value |1>')
            args[0] = nlist[i]
            args[1] = target
            time1[i] = G.test(args)
        IO.timeplotn(nlist, time1)

        # --- Time for different Fock value targets ---
        n = 10
        tarlist = [i for i in range(0,200)]
        time2 = np.zeros(len(tarlist))
        for i in range(len(tarlist)):
            print('Running Grovers for 10 qubits for Fock value |' + str(tarlist[i]) + '>')
            args[0] = n
            args[1] = tarlist[i]
            time2[i] = G.test(args)
        IO.timeplottar(tarlist, time2)

    else:
        print('\nThis is not a valid option.')
        sys.exit()

main()
