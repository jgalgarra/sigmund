import b_sim,os,sys,getopt
import numpy as np
from time import time
import bino_find_boundary

'''  Ejecuta una serie de simulaciones de bino_find_boundary y calcula la media
Ejemplo de llamada: 
bino_boundary_average.py nexper_1x1 Logistic_abs numpuntos Binary 1.0 0.1 numintentos'''

 
def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:v", ["help", "output="])
        inputfile= args[0]
        algor = args[1]
        numpoints = int(args[2])
        probtype = args[3]
        probvalue=float(args[4])
        probsd = float(args[5])
        numintentos=int(args[6])
        
    except getopt.GetoptError as err:
        print(err) # will print something like "option -a not recognized"
        print("Usage: bino_boundary_average.py nexper_1x1 Logistic_abs numpuntos Binary 1.0 0.1 numintentos")
        sys.exit(2)
    else:
        print ('Simulacion inputfile:%s algor:%s numpoints:%d probtype:%s probvalue:%f probsd:%f numintentos:%d' %\
                (inputfile,algor,numpoints,probtype,probvalue,probsd,numintentos))
        haymut=1
        haypred=0
        haysup=0
        data_save = 0
        el=0
        plants_extinction={}
        pols_extinction={} 
        #plants_extinction={'period':365,'spike':1,'start':0.05,'rate':-0.25,'numperiod':12,'species':[0,1,2]}
        #pols_extinction={'period':365,'spike':0.7,'start':0.35,'rate':-0.25,'numperiod':3,'species':[0,1]}
        dirsal='output/'
        dirent='input/'
        dirstat='output_stat_exper/'
        dirs=os.path.dirname(dirsal)
        outputdatasave=''
        osfx=''
        fichr = ''

        for i in range(numintentos):
            os.system("bino_find_boundary.py "+str(inputfile)+" "+str(algor)+" "+str(numpoints)+" "+str(probtype)+\
                                            " "+str(probvalue)+" "+str(probsd))
            probstr = str(probvalue)
            nfsal = "bd_"+str(inputfile)+"_"+str(algor)+"_"+str(probtype)+"_"+probstr+"_"+str(probsd)
            print(nfsal)
            cmd ="move "+dirstat.replace('/','\\')+str(nfsal)+".txt"+" "+dirstat.replace('/','\\')+"averages\\"+str(nfsal)+"_NUM"+str(i)+".txt"
            print(cmd)
            os.system(cmd)

        
if __name__ == "__main__":
    main()