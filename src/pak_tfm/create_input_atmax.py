#!/usr/bin/env python3
''' This script runs a simulation, takes the final populations and creates
new input matrixes with those as initial values. User must check visually
in a previous run that system reached maxima.

Conditions: Simulation file name must be xxxx_YYYYYY_zzz.sim where xxxx 
is the name of input file xxxx'''

import subprocess
com_base =  "python sigmund_standalone.py -simfile exp18_Verhulst__10.sim -y 50 -stop -rlink "
for i in range(0,no_experiments):
    #rrate = (1.0/no_experiments)*i
    rrate = 0.7
    b = subprocess.check_output(com_base+str(rrate), shell=True)
    if str(b).find("EXTINCTION") == -1:
      print("remove links "+str(rrate)+" survival")
      survival_success += 1
    else:
      print("remove links "+str(rrate)+" extinction")
print("Experiments "+str(no_experiments)+"/ Network survived " + str(survival_success)+" times")
