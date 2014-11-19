#!/usr/bin/env python3
# Example
# ./survival_count.py -simfile nar_exp1_weak_3_MAX_Verhulst__1000.sim -numexper 1 -years 200 -Bssvarsdini 20 -Bssvarsdfin 52 -Bssvarsdstep 3 -pl_species ALL
import subprocess
import argparse
parser = argparse.ArgumentParser()

parser.add_argument('-simfile', action='store', dest='simfile_name',
                    help='Simulation parameters file')

parser.add_argument('-years', action='store', dest='sim_years', default = '',
                    help='Simulation span in years')
                    
parser.add_argument('-numexper', action='store', dest='num_exper', default = '',
                    help='Number of experiments by simulation')             
                    
parser.add_argument('-Bssvarsdini', action='store', dest='Bssvarsd_ini', default = '',
                    help='Bssvarsdini')  
                    
parser.add_argument('-Bssvarsdfin', action='store', dest='Bssvarsd_fin', default = '',
                    help='Bssvarsdfin')              
                    
parser.add_argument('-Bssvarsdstep', action='store', dest='Bssvarsd_step', default = '',
                    help='Bssvarsdstep')
                    
parser.add_argument('-pl_species', action='store', dest='pl_species',
                     help='Blossom vairability, affected species')                    

conditions = parser.parse_args()


if len(conditions.simfile_name)>0:
    simfile = conditions.simfile_name
if len(conditions.sim_years)>0:
    year_periods = conditions.sim_years
if len(conditions.num_exper)>0:
    number_experiments = conditions.num_exper
if len(conditions.Bssvarsd_ini)>0:
    Bssvarsdini = conditions.Bssvarsd_ini
if len(conditions.Bssvarsd_fin)>0:
    Bssvarsdfin = conditions.Bssvarsd_fin
if len(conditions.Bssvarsd_step)>0:
    Bssvarsdstep = conditions.Bssvarsd_step
if len(conditions.pl_species)>0:
    pl = conditions.pl_species

com_base =  "python sigmund_standalone.py -simfile "+ simfile +\
            " -stop -Bssvarper 0.1 -years "+ year_periods +\
            " -Bssvartype None -Bssvarspecies "+pl

for j in range(int(Bssvarsdini),int(Bssvarsdfin),int(Bssvarsdstep)):
    survival_success = 0
    for i in range(0,int(number_experiments)):
        comando = com_base + " -Bssvarsd "+str(j/1000)
        b = subprocess.check_output(comando, shell=True)
        if str(b).find("EXTINCTION") == -1:
            survival_success += 1
#            print("SURVIVED. Survival rate "+str(survival_success/(i+1)))
#        else:
#            print("EXTINCTION. Survival rate "+str(survival_success/(i+1)))
    print(str(j/10000)+"\t"+number_experiments+"\t" + str(survival_success))