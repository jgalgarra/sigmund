'''
Created on 22/07/2015

@author: algarra
'''

import math
import numpy as np
from time import time

global f1_,f_2

def approx_binomial(n, p, size=None):
    gaussian = np.random.normal(n*p, n*p*(1-p), size=size)
    # Add the continuity correction to sample at the midpoint of each integral bin.
    gaussian += 0.5
    if size is not None:
        binomial = gaussian.astype(np.int64)
    else:
        # scalar
        binomial = int(gaussian)
    return binomial
 
def pow_mclaurin(Ranual, invperiod):
    return (1+Ranual)**invperiod -1

def calc_r_periodo_vh(Ranual, invperiod):
    return(math.pow(1 + Ranual, invperiod) - 1)

def dummy():
    calc_r_periodo_vh(0.02, invperiod)

if __name__ == '__main__':
    tini = time() 
    print("inicio")
    invperiod = 1/365
#     f_1 = 2**invperiod
#     f_2 = (2*invperiod)**(invperiod-1)
#     
#  
#     
    sumdif = 0
    Ranual = 0.02
    for i in range(50000):
        rp = calc_r_periodo_vh(Ranual, invperiod)
        #rpmc = pow_mclaurin(Ranual, invperiod)

        math.pow(1 + Ranual, invperiod) - 1
        z = -math.expm1(-rp)
        #bino = np.random.binomial(1000,z)
        apbin = approx_binomial(1000, z)
        #sumdif += abs(bino-apbin)
    #print("sumdif",sumdif)


    tfin = time()

    tfin = time()
    print("Elapsed time %.02f s" % (tfin - tini))
    