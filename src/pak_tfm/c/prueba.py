from helloworld import helloworld
from time import time
import math

lim = 10000000

j = 0
tinic=time()
while (j<lim):
    r = helloworld(1.0)
    j+=1
tfin=time()
print("Con C. Elapsed time %.02f s" % (tfin-tinic))

j=0
tinic=time()
while (j<lim):
    r = math.expm1(1.0)
    j+=1
tfin=time()
print("Python solo. Elapsed time %.02f s" % (tfin-tinic))