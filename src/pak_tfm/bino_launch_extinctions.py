import b_sim,os,sys,getopt,time
import numpy as np

''' Ejemplo de llamada: 
bino_launch_extinctions.py nexper_1x1 80 100 2 1 Logistic_abs Binary 0.01 '''

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
        pbinic = float(args[1])
        pbfin = float(args[2])
        numpasos = int(args[3])
        numsims = int(args[4])
        algor = args[5]
        btype = args[6]
        pb_sd = float(args[7])
        
    except getopt.GetoptError as err:
        print(err) # will print something like "option -a not recognized"
        print("Usage: bino_launch_extinctions.py inputfile pbinic pbfin numpasos numsims algor btype pb_sd")
        sys.exit(2)
    else:
        print ('Simulacion inputfile:%s pbinic:%f pbfin:%f numpasos:%d numsims:%d  algor:%s  btype:%s pb_sd:%f' % (inputfile,pbinic,pbfin,numpasos,numsims,algor,btype,pb_sd))
        displayinic=0
        ciclos=100
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
        dirs=os.path.dirname(dirsal)
        outputdatasave=''
        osfx=''
        fichr = ''
        
        incpb=(pbfin-pbinic)/numpasos
        pbiter =pbinic
        z = []
        paso = 1
        while (pbiter <= pbfin):
            print("Step "+str(paso)+". pbiter "+str(pbiter)+" "+time.asctime( time.localtime(time.time()) ))
            sumaextinct = 0
            pb = pbiter/100
            for j in range(0,numsims):
                Na,Nb,Nc,Ra,Rb,RequA,RequB,maxa,maxb,maxreff,minreff,maxequs,minequs,extinction =b_sim.bino_mutual(inputfile,ciclos,haypred,haysup,outputdatasave,dirs,\
                                                                     dirent,dirsal,eliminarenlaces=el,pl_ext=plants_extinction,pol_ext=pols_extinction,\
                                                                     os=osfx,fichreport=fichr,com='',algorithm=algor,plants_blossom_prob=pb,\
                                                                     plants_blossom_sd=pb_sd,plants_blossom_type=btype,blossom_pert_list=['ALL'],\
                                                                     verbose=False,exit_on_extinction=True)        
                sumaextinct+=extinction
            z.append([pb,sumaextinct/numsims])
            paso+=1
            print("ZETA "+str(z))
            pbiter += incpb
        print("ZETA FINAL "+str(z))
        if (btype == 'Binary'):
            nfescritura = 'output_stat_exper/z_'+inputfile+'_'+btype+'_'+algor+'_ini'+str(pbinic)+'_fin'+str(pbfin)+'_numsims'+str(numsims)+'_pasos_'+str(numpasos)+'.txt'
        else:
            nfescritura = 'output_stat_exper/z_'+inputfile+'_'+btype+algor+'_sd'+str(pb_sd)+'_ini'+str(pbinic)+'_fin'+str(pbfin)+'_numsims'+str(numsims)+'_pasos_'+str(numpasos)+'.txt'  
        dlmwritesim(nfescritura,z)
  
        return

if __name__ == "__main__":
    main()