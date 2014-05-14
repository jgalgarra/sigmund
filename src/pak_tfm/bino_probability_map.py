import b_sim,os,sys,getopt
import numpy as np
from time import time
import bino_probability_slice

''' Ejemplo de llamada: 
bino_probability_map.py nexper_1x1 fileboundary 10 1 100 Logistic_abs '''

def dlmwritesim(inputname,Nin):
    print ("Output file %s" % inputname)
    salida=open(inputname,'w',encoding='utf-8')
    for linea in Nin:
        for i in range(len(linea)-1):
            y="%.08f" % linea[i]
            salida.write(y+',')
        salida.write(str(linea[-1])+'\n');
    salida.close()
    return


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:v", ["help", "output="])
        inputfile= args[0]
        fileboundary = args[1]
        numpoints = int(args[2])
        numexperiments = int(args[3])
        ciclos=int(args[4])
        algor = args[5]
        
    except getopt.GetoptError as err:
        print(err) # will print something like "option -a not recognized"
        print("Usage: bino_probability_map.py nexper_1x1 fileboundarynumpoints numexperiments numyears algor")
        sys.exit(2)
    else:
        print ('Simulacion inputfile:%s boundary file:%s numexperiments:%d numyears:%d algor:%s direction:%s' % (inputfile,fileboundary,numpoints,numexperiments,ciclos,algor))
       
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

        btype = 'Binary'
        pb_sd = 0.01
        pb = 1.00

        def filerow_to_point(lstent):
            return int(lstent[0].split(',')[0]),int(lstent[0].split(',')[1])

        l_boundary=b_sim.dlmreadlike(fileboundary,dirstat)
        print(l_boundary)
        numpuntos = len(l_boundary)
        pcentral = l_boundary[numpuntos//2]
        pcx,pcy = filerow_to_point(pcentral)
        for puntolista in l_boundary:
            puntox, puntoy = filerow_to_point(puntolista)
            print(puntox,puntoy)
            if (puntox<pcx):
                sentido = 'H'
            elif (puntox>pcx):
                sentido = 'V'
            os.system("bino_probability_slice.py "+str(inputfile)+" "+str(puntox)+" "+str(puntoy)+" "+str(numpoints)+\
                                            " "+str(numexperiments)+" "+str(ciclos)+" "+str(algor)+" "+str(sentido))
            if (puntox == pcx):
                os.system("bino_probability_slice.py "+str(inputfile)+" "+str(puntox)+" "+str(puntoy)+" "+str(numpoints)+\
                                            " "+str(numexperiments)+" "+str(ciclos)+" "+str(algor)+" H")
                os.system("bino_probability_slice.py "+str(inputfile)+" "+str(puntox)+" "+str(puntoy)+" "+str(numpoints)+\
                                            " "+str(numexperiments)+" "+str(ciclos)+" "+str(algor)+" V")
            
        return
        '''
        if (direction=='H'):
            vcoord = np.linspace(xcentral-numpoints,xcentral+numpoints,numpoints*2+1)
            print(vcoord)
            indexvar = 0
        else:
            vcoord = np.linspace(ycentral-numpoints,ycentral+numpoints,numpoints*2+1)
            print(vcoord)
            indexvar = 1
        listaresults = []
        tinic_script=time()
        firstnonzerofound = False
        for i in (vcoord):
            tinic=time()
            cuentasurvival = 0
            print("Point:",int(i),y)
            if (len(listaresults)>=2) and (listaresults[-1][2]+listaresults[-2][2]==2.0):
                listaresults.append([int(i),y,1.0]) if direction=='H' else listaresults.append([x,int(i),1.0])
            else:   
                if not(firstnonzerofound):
                    exper = min(100,round(numexperiments/2))
                else:
                    exper = numexperiments
                if direction=='H':
                    N0plants_val = str(int(i))
                    N0pols_val=str(y)  
                else:
                    N0plants_val=str(x)
                    N0pols_val=str(int(i))
                for j in range(exper):
                    Na,Nb,Nc,Ra,Rb,RequA,RequB,maxa,maxb,maxreff,minreff,maxequs,minequs,extinction =b_sim.bino_mutual(inputfile,ciclos,haypred,haysup,outputdatasave,dirs,\
                                                                         dirent,dirsal,eliminarenlaces=el,pl_ext=plants_extinction,pol_ext=pols_extinction,\
                                                                         os=osfx,fichreport=fichr,com='',algorithm=algor,plants_blossom_prob=pb,\
                                                                         plants_blossom_sd=pb_sd,plants_blossom_type=btype,blossom_pert_list=['ALL'],\
                                                                         verbose=False,exit_on_extinction=True,N0plants=N0plants_val,N0pols=N0pols_val)        
                    cuentasurvival += (extinction == False)
                listaresults.append([int(i),y,float(cuentasurvival/exper)]) if direction=='H' else listaresults.append([x,int(i),float(cuentasurvival/exper)]) 
                if not(firstnonzerofound):
                    firstnonzerofound = (cuentasurvival>0)
            tfin = time()
            print ("Elapsed time %.02f s" % (tfin-tinic))
            print(listaresults)
            print()
        tfin_script = time()
        print ("Total time of simulation %.02f s" % (tfin_script-tinic_script))
           
        nfescritura = 'output_stat_exper/probslice_'+inputfile+'_'+btype+'_'+algor+'_xcentral'+str(xcentral)+'_ycentral'+str(ycentral)+'_numpoints'+str(numpoints)+'_nexper_'+str(numexperiments)+'_years_'+str(ciclos)+'_'+str(direction)+'.txt'
        dlmwritesim(nfescritura,listaresults)
        return'''

if __name__ == "__main__":
    main()