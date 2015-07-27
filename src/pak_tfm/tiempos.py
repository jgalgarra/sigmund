'''
Created on 22/07/2015

@author: algarra
'''

import math
import numpy as np
from time import time

global f1_,f_2
 
def pow_mclaurin(Ranual, invperiod):
    return (1+Ranual)**invperiod -1

def calc_r_periodo_vh(Ranual, invperiod):
    return(math.pow(1 + Ranual, invperiod) - 1)

def dummy():
    calc_r_periodo_vh(0.02, invperiod)

if __name__ == '__main__':
    tini = time() 
    print("inicio")
#     invperiod = 1/365
#     f_1 = 2**invperiod
#     f_2 = (2*invperiod)**(invperiod-1)
#     
#  
#     
#     sumdif = 0
#     Ranual = 0.02
#     for i in range(500000):
#         rp = calc_r_periodo_vh(Ranual, invperiod)
#         rpmc = pow_mclaurin(Ranual, invperiod)
#         sumdif += abs(rp-rpmc)
#         #math.pow(1 + Ranual, invperiod) - 1
#         z = -math.expm1(-rp)
#         np.random.binomial(1000,z)
#         #np.random.binomial(1000, 0.0001)
#     print("sumdif",sumdif)

    
    [abs(-3.09876) for i in range(10000000)]
    tfin = time()

    tfin = time()
    print("Elapsed time %.02f s" % (tfin - tini))
    