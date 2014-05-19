import b_sim
import os,sys,getopt,time
import numpy as np
import math

global num
global probtype
global probvalue
global probsd

''' Ejemplo de llamada: 
bino_find_boundary.py nexper_1x1 Logistic_abs numpuntos Binary 1.0 0.1'''

def dlmwritesim(inputname,Nin):
    print ("Output file %s" % inputname)
    salida=open(inputname,'w',encoding='utf-8')
    for linea in Nin:
        for i in range(len(linea)-1):
            y="%d" % linea[i]
            salida.write(y+',')
        salida.write(str(linea[-1])+'\n');
    salida.close()
    return



def search_transition_candidate(inputfile, Nplants, Npols, algor, ciclos, haypred, haysup, dirsal, dirent, dirs, outputdatasave,fichr, test_punto, extinction, paso):
    cuentasurv = 0
    limsurv = 6 if (paso>1) else 20
    while cuentasurv <= limsurv:
        Na, Nb, Nc, Ra, Rb, RequA, RequB, maxa, maxb, maxreff, minreff, maxequs, minequs, extinction = \
        b_sim.bino_mutual(inputfile, ciclos, haypred, haysup, outputdatasave, dirs, 
            dirent, dirsal, 
            fichreport='', com='', algorithm=algor, plants_blossom_prob=probvalue,\
                plants_blossom_sd= probsd, plants_blossom_type = probtype, blossom_pert_list=['ALL'], 
            verbose=False, exit_on_extinction=True, N0plants=Nplants, N0pols=Npols)
        if extinction:
            cuentasurv = 0
            #print ("System Collapse", test_punto)
            return(True)                             #returns extinction
        else:
            cuentasurv = cuentasurv + 1
    return(False)


#def find_transition_candidate(inputfile, N0plants_arg, N0pols_arg, algor, ciclos, haypred, haysup, dirsal, dirent, dirs, outputdatasave, fichr, test_punto, paso):
def find_transition_candidate(inputfile, algor, ciclos, haypred, haysup, dirsal, dirent, dirs, outputdatasave, fichr, test_punto, paso, index_interval):
 
    extinction = True
    while extinction:
        extinction = search_transition_candidate(inputfile, str(int(test_punto[0])), str(int(test_punto[1])), algor, ciclos, haypred, haysup, dirsal, dirent, dirs, outputdatasave, fichr, test_punto, extinction, paso)
        if (extinction):
            test_punto[index_interval] = round(test_punto[index_interval] + paso)
        else:
            print("System Survival", test_punto)
    return(test_punto)


def find_transition_point(inputfile, algor, ciclos, haypred, haysup, dirsal, dirent, dirs, outputdatasave, fichr, rango_N_ini, rango_N_fin, constant_coord, interval_direction):
    if (interval_direction=='Horizontal'):
        test_punto = [rango_N_ini, constant_coord]
        index_interval = 0
    else:
        test_punto = [constant_coord,rango_N_ini]
        index_interval = 1
    paso = max([1 , round(abs(rango_N_ini - rango_N_fin) / 20)]);
    punto_found = False
    while not punto_found:
        print ("Paso:", paso)
        punto_transicion = find_transition_candidate(inputfile, algor, ciclos, haypred, haysup, dirsal, dirent, dirs, outputdatasave, fichr, test_punto, paso, index_interval)
        print ("Candidato a punto transicion:", punto_transicion)
        if (paso == 1):
            punto_found = True
        else:
            rango_N_ini = max(1,round(test_punto[index_interval] - 2 * paso))
            rango_N_fin = test_punto[index_interval]
            paso = max(1, round(abs(rango_N_ini - rango_N_fin) / 20))
            if (interval_direction=='Horizontal'):
                test_punto = [rango_N_ini, constant_coord] 
            elif (interval_direction=='Vertical'):  
                test_punto = [constant_coord, rango_N_ini]
    return test_punto

def calc_val_verhuslt(r1,r2,c1,c2,alpha1,alpha2,b12,b21):
    k1_A = (c2*b21*alpha1+c1*b12*b21)
    k1_B = (alpha1*alpha2 + c1*b12*r2 - c2*b21*r1 -b12*b21)
    k1_C = -1*(r1*alpha2 + b12* r2)

    raiz_k1_1 = (-1*k1_B+math.sqrt(k1_B**2-4*k1_A*k1_C))/(2*k1_A)
    raiz_k1_2 = (-1*k1_B-math.sqrt(k1_B**2-4*k1_A*k1_C))/(2*k1_A)

    k2_A = (c1*b12*alpha2+c2*b21*b12)
    k2_B = (alpha2*alpha1 + c2*b21*r1 - c1*b12*r2 -b21*b12)
    k2_C = -1*(r2*alpha1 + b21* r1)

    raiz_k2_1 = (-1*k2_B+math.sqrt(k2_B**2-4*k2_A*k2_C))/(2*k2_A)
    raiz_k2_2 = (-1*k2_B-math.sqrt(k2_B**2-4*k2_A*k2_C))/(2*k2_A)
    ssup = [raiz_k1_1,raiz_k2_1]
    sinf = [raiz_k1_2,raiz_k2_2]
    return sinf,ssup

def main():
    global probtype,probvalue, probsd 
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:v", ["help", "output="])
        inputfile= args[0]
        algor = args[1]
        numpuntos = int(args[2])
        probtype = str(args[3])
        probvalue = float(args[4])
        probsd = float(args[5])
        
    except getopt.GetoptError as err:
        print(err) # will print something like "option -a not recognized"
        print("Usage: bino_find_boundary.py inputfile algor 10 Binary 1.0 0.1")
        sys.exit(2)
    else:
        print ('Simulacion inputfile:%s algor:%s probtype:%s probvalue:%.02f probsd:%.02f' % (inputfile,algor,probtype,probvalue,probsd))
        ciclos=100
        haypred=0
        haysup=0
        el=0
        plants_extinction={}
        pols_extinction={} 
        dirsal='output/'
        dirent='input/'
        dirs=os.path.dirname(dirsal)
        outputdatasave=''
        fichr = ''
        filename_a=inputfile+'_a.txt'
        l_minputchar_a=b_sim.dlmreadlike(filename_a,'input/')
        rplant0 = float(l_minputchar_a[-2][0])-float(l_minputchar_a[-1][0])
        #Nplant0_conf = int(l_minputchar_a[-4][0])
        bplant0 = float(l_minputchar_a[0][0])
        filename_b=inputfile+'_b.txt'
        l_minputchar_b=b_sim.dlmreadlike(filename_b,'input/')
        rpol0 = float(l_minputchar_b[-2][0])-float(l_minputchar_b[-1][0])
        bpol0 = float(l_minputchar_b[0][0])
        #Npol0_conf = int(l_minputchar_b[-4][0])
        if (algor!='Verhulst'):
            Kplant0 = int(l_minputchar_a[-3][0])
            Ncritpol = int(math.ceil(-rplant0/bplant0))
            Kpol0 = int(l_minputchar_b[-3][0])
            Ncritplant = int(math.ceil(-rpol0/bpol0))
        else:
            alpha_plant0 = float(l_minputchar_a[-3][0])
            c_alpha_plant0 = float(l_minputchar_a[-4][0])
            alpha_pol0 = float(l_minputchar_b[-3][0])
            c_alpha_pol0 = float(l_minputchar_b[-4][0])
            bplant0 = float(l_minputchar_a[0][0])
            sinf,ssup = calc_val_verhuslt(rplant0,rpol0,c_alpha_plant0,c_alpha_pol0,alpha_plant0,alpha_pol0,bplant0,bpol0)
            #print(sinf,ssup)
            Kplant0 = math.ceil(ssup[0])
            Kpol0 = math.ceil(ssup[1])
            Ncritplant = math.ceil(sinf[0])
            Ncritpol = math.ceil(sinf[1])
        
        print("Plant r:%.02f b:%.06f NcritPlant:%.d K:%d" % (rplant0,bplant0,Ncritplant,Kplant0))
        print("Pollinator r:%.02f b:%.06f NcritPol:%.d K:%d" % (rpol0,bpol0,Ncritpol,Kpol0))

        ''' Finding corner point '''
        if (probtype == 'Binary') and (probvalue == 1.0):
            cornerpoint = [Ncritplant,Ncritpol]
        else:
            test_punto = find_transition_point(inputfile, algor, ciclos, haypred, haysup,\
                                               dirsal, dirent, dirs,\
                                               outputdatasave, fichr, Ncritpol-20, Ncritpol+20,\
                                               Ncritplant, interval_direction = 'Vertical') 
            cornerpoint = [int(test_punto[0]),int(test_punto[1])]
            
        print("Corner point:%d,%d" % (cornerpoint[0],cornerpoint[1]))
        
        numpuntos = numpuntos//2
        puntospolcalculo = np.round(np.linspace(Kpol0,cornerpoint[1],numpuntos+1))
        print(puntospolcalculo)
        
        puntosgraf = [] 
        repeat_count = 0
        margen = 0
        rango_N_ini = round(cornerpoint[0]/4)
        rango_N_fin = cornerpoint[0]
        for i in puntospolcalculo[0:-1]:
            repetir = True
            while (repetir):
                if repeat_count >=3:
                    print("ALARM!!. Exit")
                    exit()
                test_punto = find_transition_point(inputfile, algor, ciclos, haypred, haysup,\
                                                   dirsal, dirent, dirs,\
                                                   outputdatasave, fichr, rango_N_ini-margen, rango_N_fin,\
                                                   int(i), interval_direction = 'Horizontal')
                if (len(puntosgraf)>1) and (puntosgraf[-1][0]>=test_punto[0]):
                    repetir = True
                    print("WARNING!. Repetition !!!")
                    repeat_count+=1
                    margen = rango_N_ini//2
                else:
                    puntosgraf.append([int(test_punto[0]),int(test_punto[1])])
                    print(puntosgraf)
                    rango_N_ini = puntosgraf[-1][0]
                    repetir = False
                    repeat_count =0
                    margen = 0
        '''if (len(puntosgraf)>1):            
            test_punto = find_transition_point(inputfile, algor, ciclos, haypred, haysup,\
                                                   dirsal, dirent, dirs, outputdatasave, fichr,\
                                                   cornerpoint[1], puntosgraf[-1][1], round(cornerpoint[0]-puntosgraf[-1][0])/2, interval_direction = 'Vertical')'''
                    
        puntosgraf.append(cornerpoint)
        
        puntosplcalculo = np.round(np.linspace(cornerpoint[0],Kplant0,numpuntos+1))
        print(puntosplcalculo)
        rango_N_ini = max(1,round(cornerpoint[1]//4))
        rango_N_fin = cornerpoint[1]
        margen = 0
        for i in puntosplcalculo[1:]:
            repeat_count = 0
            repetir = True
            while (repetir):
                if repeat_count >=2:
                    print("ALARM!!. Exit")
                    exit()
                test_punto = find_transition_point(inputfile, algor, ciclos, haypred, haysup,\
                                                   dirsal, dirent, dirs, outputdatasave, fichr,\
                                                   rango_N_ini-margen, rango_N_fin, int(i), interval_direction = 'Vertical')
                if (len(puntosgraf)>1) and (puntosgraf[-1][1]<=test_punto[1]):
                    repetir = True
                    print("WARNING!. Repetition !!!")
                    repeat_count +=1
                    margen = rango_N_ini//2
                else:
                    puntosgraf.append([int(test_punto[0]),int(test_punto[1])])
                    print(puntosgraf)
                    rango_N_ini = int(puntosgraf[-1][1]/2)
                    repetir = False
                    repeat_count =0
                    margen = 0
            
        #tailstr = "_Pl_r_b_%0.6f_%0.6f_Pol_r_b_%.06f_%0.6f" % (rplant0,bplant0,rpol0,bpol0)
        nfescritura = 'output_stat_exper/bd_'+inputfile+'_'+algor+'_'+probtype+'_'+str(probvalue)+'_'+str(probsd)+'.txt'  
        dlmwritesim(nfescritura,puntosgraf)
        
        return

if __name__ == "__main__":
    main()