'''
Created on 10/05/2012

@author: Javier Garcia Algarra

This module includes the mutualist algorithm
'''
import numpy as np
from time import time
import math
import datetime
#import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.backends.backend_tkagg
import os
import cProfile

global ancho
global alto
global resolucion
global pathout
global hayreport
global frep
global output_verbose
global invperiod
global cuentaperpl,cuentaperpol
global Logistic_abs
global model_r_alpha
global tol_extincion
global hay_bssvar
global pendiente,sd,periodo
global diasdelanio
global ax

diasdelanio = 365

ancho=16
alto=10
resolucion=600
invperiod = 1/diasdelanio
tol_extincion = -0.0001

'''signfunc = lambda x : (x>0) - (x<0)'''

def signfunc(x):
    '''if (x>=0):
        return(1)
    else:
        return(-1)'''
    '''return(math.copysign(1, x))'''
    if abs(x)==x:
        return(1)
    else:
        return(-1)

def calc_coef_May(minputchar,r,K):
    return(minputchar if r!=0 else 0)

def calc_r_periodo(Ranual,invperiod):
    #return(math.expm1(Ranual*invperiod))   # Since * is faster than div we multiply by invperiod = 1/period
    return(math.pow(1+Ranual,invperiod) - 1)

def calc_r_periodo_vh(Ranual,invperiod):
    #return(math.expm1(Ranual/period)) 
    return(math.pow(1+Ranual,invperiod) - 1)

def cuentaenlaces(mat_in):
    cuenta=0;
    fil=len(mat_in)-3
    col=len(mat_in[0])
    for i in range(fil):
        for j in range(col):
            if abs(mat_in[i][j])>0:
                cuenta+=1
    return(cuenta)

def borraenlace(mat_in):
    intentos=0
    fil=len(mat_in)-5
    col=len(mat_in[0])
    rfil = (np.random.random_integers(0,fil-1))
    rcol = (np.random.random_integers(0,col-1))
    while (mat_in[rfil][rcol]==0) and (intentos<fil*col):
        rfil = (np.random.random_integers(0,fil-1))
        rcol = (np.random.random_integers(0,col-1))
        intentos+=1
    mat_in[rfil][rcol]=0
    return(mat_in,rfil,rcol)

def dlmreadlike(inputfile,direntrada):
    try:
        data = open(direntrada+inputfile)
        mtx_input=[]
        for each_line in data:
            try:
                #print (each_line.split('\t'))
                mtx_input.append(each_line.strip().split('\t'))
            except ValueError:
                pass
        data.close()
        return(mtx_input)
    except IOError:
        print('The datafile %s is missing!' %inputfile)
        return(0)

def dlmwritelike(inputfile,nperiod,Nin,dirsalida,os):
    dsal = dirsalida.replace('\\','/')
    nsal= 'output_data_'+inputfile+'_'+os+'_'+str(nperiod)+'.txt'
    print ("Output file %s" % dsal+nsal)
    salida=open(dsal+nsal,'w',encoding='utf-8')
    for linea in Nin:
        for i in range(len(linea)-1):
            y="%.012f" % linea[i]
            #salida.write(str(linea[i])+',')
            salida.write(y+',')
        salida.write(str(linea[-1])+'\n');
    salida.close()
    return(nsal)

def informa(text):
    global hayreport
    global frep
    global output_verbose
    
    if output_verbose:
        '''print(text,file=sys.stdout)'''
        print(text)
        pass
    if hayreport: 
        frep.write(text+'<br>')
    return(0)

def val_mutMay(r_species,beta,period,N1,N2,K1):
    rspneq = calc_r_periodo_vh(abs(r_species),invperiod)
    betaMay = beta*K1/abs(r_species)
    termEq = round(betaMay*N1*N2/K1)
    rMay = betaMay*N2/K1
    if (abs(termEq)> 1):
        #incEq = np.random.binomial(round(abs(termEq)),1-math.exp(-rspneq))
        incEq = np.random.binomial(round(abs(termEq)),-math.expm1(-rspneq)) 

    else:
        incEq = 0
    ret = [incEq,rMay]
    return(ret)

def ciclo_May(r_species, rM, period, inctermMay, Nindivs, K):
    rspneq = calc_r_periodo_vh(abs(r_species),invperiod)
    signosp = signfunc(r_species)
    termEq = Nindivs * (signosp-(Nindivs/K))
    rcal = r_species*((1-(Nindivs/K)) + rM)
    if (abs(termEq)> 1):
        #incEq = np.random.binomial(round(abs(termEq)),1-math.exp(-rspneq))
        incEq = np.random.binomial(round(abs(termEq)),-math.expm1(-rspneq))
    else:
        incEq = 0
    ret = [Nindivs + signfunc(termEq)*incEq + abs(inctermMay),signfunc(termEq)*rcal]
    return(ret)
    #return(Nindivs + signfunc(term3Eq)*inc3Eq)
               
def ciclo_verhulst(rtot_species, reqsum, period, Nindivs, cAlpha, Alpha):
    roz = (Alpha+cAlpha*reqsum)*Nindivs
    rcal = rtot_species - roz
    #print(rcal,reqsum,roz,Alpha)
    '''if (rtot_species>0):
        rcal -= (rtot_species * Nindivs*invK) 
    elif L_abs :       # |reff| correction
        rcal += (rtot_species*Nindivs*invK)'''      
    rspneq = calc_r_periodo_vh(rcal,invperiod) if (rcal>=0) else calc_r_periodo_vh(-rcal,invperiod)
    incNmalth= np.random.binomial(Nindivs,-math.expm1(-rspneq))
    if (rcal>0):
        inc_pop = Nindivs + incNmalth 
    else:
        inc_pop = Nindivs - incNmalth
    return([inc_pop,rcal])     


def ciclo_new_model(rtot_species, period, Nindivs, invK, L_abs):
    rcal = rtot_species
    if (rtot_species>0):
        rcal -= (rtot_species * Nindivs*invK) 
    elif L_abs :       # |reff| correction
        rcal += (rtot_species*Nindivs*invK)      
    rspneq = calc_r_periodo_vh(rcal,invperiod) if (rcal>=0) else calc_r_periodo_vh(-rcal,invperiod)
    incNmalth= np.random.binomial(Nindivs,-math.expm1(-rspneq))
    if (rcal>0):
        inc_pop = Nindivs + incNmalth 
    else:
        inc_pop = Nindivs - incNmalth
    return([inc_pop,rcal]) 

def ciclo_none(rtot_species, period, Nindivs, invK, L_abs):
    rcal = rtot_species
    rspneq = calc_r_periodo_vh(rcal,invperiod) if (rcal>=0) else calc_r_periodo_vh(-rcal,invperiod)
    incNmalth= np.random.binomial(Nindivs,-math.expm1(-rspneq))
    if (rcal>0):
        inc_pop = Nindivs + incNmalth 
    else:
        inc_pop = Nindivs - incNmalth
    return([inc_pop,rcal])    

def calc_mutualism_params(minputchar_p,Alpha_p,Nindividuals_p,Nindividuals_q,r_eqsum,term_May,rMay,k,j,n,r_p,r_muerte,days_year):
    coef_bij_matrix=float(minputchar_p[j][n]*Nindividuals_q[k][j])
    lr_eqsum = float(r_eqsum + coef_bij_matrix)
    retMay = val_mutMay(r_p[n]-r_muerte,calc_coef_May(minputchar_p[j][n],r_p[n]-r_muerte,Alpha_p[n]),days_year,Nindividuals_p[k][n],Nindividuals_q[k][j],Alpha_p[n])
    lterm_May = term_May+retMay[0]
    lrMay = rMay + retMay[1]
    return lr_eqsum,lterm_May,lrMay


def init_lists_pop(periods,numspecies_p):
    Alpha_p = []
    cAlpha_p = []
    rowNindividuals_p = []
    r_p = []
    rd_p = []
    Nindividuals_p=np.zeros([periods,numspecies_p])
    rp_eff = np.zeros([periods,numspecies_p],dtype=float)
    rp_equs = np.zeros([periods,numspecies_p],dtype=float)
    return rowNindividuals_p, Alpha_p, cAlpha_p, r_p, rd_p, Nindividuals_p, rp_eff, rp_equs #, maxima_ind_p, minima_ind_p


def init_params_population(r_p, rd_p, n):
    r_eqsum = term_May = rMay = p_devorados = rtot_p = 0
    #signo = signfunc(r_p[n])
    r_muerte = rd_p[n]
    pop_ini = 0
    return r_muerte, r_eqsum, term_May, rMay, rtot_p, p_devorados, pop_ini


def perturbation(pl_ext, n, rd_p, cuentaperp, inicioextp, periodoextp, spikep, k):
    r_m = rd_p
    cuentaper = cuentaperp
    if (k >= inicioextp) and (n in pl_ext['species']):
        posindextp = (k - inicioextp) % periodoextp
        if (posindextp >= 0) and (posindextp < spikep):
            r_m = rd_p + pl_ext['rate']
        if (posindextp == spikep - 1) and (n==pl_ext['species'][0]):
            cuentaper = cuentaper + 1
    return r_m, cuentaper


def init_perturbations(pl_ext, pol_ext, yearperiods, inicioextplantas, inicioextpolin, hayextplantas, hayextpolin):
    if hayextplantas:
        inicioextplantas = round(yearperiods * pl_ext['start'] * diasdelanio)
        nperpl = pl_ext['numperiod']
        periodoextpl = pl_ext['period']
        spikepl = round(periodoextpl * pl_ext['spike'])
        informa("Perturbations. Plants species %s, period (years): %d, numperiods: %d, spike (fraction of period): %0.2f, rate: %.03f, start (year): %.02f" % (pl_ext['species'], periodoextpl / diasdelanio, nperpl, pl_ext['spike'], float(pl_ext['rate']), pl_ext['start']))
    else:
        inicioextplantas = nperpl = periodoextpl = spikepl = 0
    if hayextpolin:
        inicioextpolin = round(yearperiods * pol_ext['start'] * diasdelanio)
        nperpol = pol_ext['numperiod']
        periodoextpol = pol_ext['period']
        spikepol = round(periodoextpol * pol_ext['spike'])
        informa("Perturbations. Pollinators species %s, period (years): %d, numperiods: %d, spike (fraction of period): %0.2f, rate: %.03f, start (year): %.02f" % (pol_ext['species'], periodoextpol / diasdelanio, nperpol, pol_ext['spike'], float(pol_ext['rate']), pol_ext['start']))
    else:
        inicioextpolin = nperpol = periodoextpol = spikepol = 0
    return nperpl, inicioextplantas, periodoextpl, spikepl, nperpol, inicioextpolin, periodoextpol, spikepol

def sinusoidal(val):
    global periodo,sd
    modifier = 0.5*(1.0001+np.sin((2*np.pi/periodo)*val))
    return np.random.normal(1,sd*modifier)

def lineal(val):
    global numanyos,pendiente
    modifier = 1+pendiente*(val/numanyos)
    return np.random.normal(1,sd*modifier)

def sinmodulacion(val):
    global sd
    return(np.random.normal(1,sd))

def calcbssvarespecie(valblossomperiod,valsd,funct,nys):
    global numanyos,sd,blossomperiod
    blossomperiod = valblossomperiod
    sd = valsd
    numanyos = nys
    indiffs = []
    valdiffs=[]
    for m in range (0,numanyos):
        pospl = funct(m)
        diff = 1-((2/blossomperiod)*(abs(1-pospl)))
        if (diff < 0):
            diff = 0
        valdiffs.append(diff)
        indiffs.append(m)     
    return indiffs,valdiffs

def predators_effect(p_devorados, days_year, j, Nindividuals_p, Nindividuals_c, minputchar_c, numspecies_c, n, k):
    for j in range(numspecies_c):
        if (minputchar_c[n][j] > 0):
            rceff = calc_r_periodo_vh(Nindividuals_c[k][j] * minputchar_c[n][j], invperiod)
            p_devorados = p_devorados + np.random.binomial(Nindividuals_p[k][n], -math.expm1(-1 * rceff))
    return p_devorados, j, rceff


def populations_evolution(n,strtype, numspecies_p, algorithm, hay_foodweb, p_ext, May, haymut,\
                          days_year, rp_eff, rp_eq, minputchar_p, j, cAlpha_p, Alpha_p, r_p, rd_p, \
                          Nindividuals_p, numspecies_q, Nindividuals_q, Nindividuals_c, minputchar_c, numspecies_c, \
                          inicioext, hayext, nper, periodoext, spike, k, model_r_a):
    ''' for n in range (numspecies_p):'''
    global cuentaperpl,cuentaperpol
    #rceff = retl = rcalc = 0
    r_muerte, r_eqsum, term_May, rMay, rtot_p, p_devorados, pop_p = init_params_population(r_p, rd_p, n)
    '''try:'''
    #if (Nindividuals_p[k][n] > 0):
    if (Nindividuals_p[k,n]):           # Much faster than (Nindividuals_p[k][n] > 0)
        # Extinciones de plantas
        if hayext:
            if (strtype=='Plant'):
                cuentaperext = cuentaperpl  
            else: 
                cuentaperext = cuentaperpol
            if (cuentaperext < nper):
                r_muerte, cuentaperext = perturbation(p_ext, n, rd_p[n], cuentaperext, inicioext, periodoext, spike, k)
                if (strtype=='Plant'):
                    cuentaperpl = cuentaperext
                else:
                    cuentaperpol = cuentaperext
        if haymut:
            r_eqsum = np.dot(minputchar_p[0:numspecies_q,n],Nindividuals_q[k,:])
        elif (May):
            for j in range(numspecies_q):
                r_eqsum, term_May, rMay = calc_mutualism_params(minputchar_p, Alpha_p, Nindividuals_p, Nindividuals_q, r_eqsum, term_May, rMay, k, j, n, r_p, r_muerte, days_year)
        rtot_p = r_p[n] + r_eqsum - r_muerte

        # Efecto de los depredadores
        if hay_foodweb:
            p_devorados, j, rceff = predators_effect(p_devorados, days_year, j, Nindividuals_p, Nindividuals_c, minputchar_c, numspecies_c, n, k)
        # New algorithm
        if (model_r_a):
            retl = ciclo_verhulst(rtot_p, r_eqsum, days_year, Nindividuals_p[k,n], cAlpha_p[n], Alpha_p[n])
        elif not (May):
            retl = ciclo_new_model(rtot_p, days_year, Nindividuals_p[k,n], 1/Alpha_p[n], Logistic_abs)
        else:
            retl = ciclo_May(r_p[n] - r_muerte, rMay, days_year, term_May, Nindividuals_p[k,n], Alpha_p[n])           
        pop_p = retl[0] - p_devorados
        if not(pop_p):
            informa("Day %d (year %d). %s species %d extincted" % (k, k//diasdelanio, strtype, n))
    Nindividuals_p[k+1][n] = pop_p
    if (pop_p):
        rp_eff[k+1][n] = retl[1]
        rp_eq[k+1][n] = rtot_p 
    else:
        rp_eff[k+1][n] =  rp_eff[k][n]
        rp_eq[k+1][n] = rp_eq[k][n]
    return 0
    
'''def calc_compatib_plantas(numspecies,probcoinc,blossom_pert_list):
    lcomp = np.ones(numspecies,dtype=float)
    if (probcoinc == 1.0):
        return lcomp
    for i in blossom_pert_list:
        lcomp[i] = np.random.binomial(n=1, p=probcoinc)
    return lcomp'''   

def calc_compatib_plantas(numspecies,probcoinc,blossom_pert_list):
    lcomp = []
    for i in range(0,numspecies):
        if i in blossom_pert_list:
            lcomp.append(np.random.binomial(n=1, p=probcoinc))
        else:
            lcomp.append(1.0)
    return np.array(lcomp)       

def calc_blossom_effect(numspecies_a,nrows_a,ncols_a,nrows_b,ncols_b,numspecies_b,plants_blossom_type,\
                        plants_blossom_prob,plants_blossom_sd,blossom_pert_list):
    if (plants_blossom_type == 'Binary'):
        lcompatibplantas = calc_compatib_plantas(numspecies_a,plants_blossom_prob,blossom_pert_list)
    else:
        lcompatibplantas = []
        for g in range(0,numspecies_a):
            if g in blossom_pert_list:
                lcompatibplantas.append(abs(np.random.normal(plants_blossom_prob, plants_blossom_sd, 1)))
            else:
                lcompatibplantas.append(1)
    minputchar_a_mask = np.ones([nrows_a,ncols_a])
    for i in range(0,numspecies_b):
        minputchar_a_mask[i,:]=lcompatibplantas
    minputchar_b_mask = np.ones([nrows_b,ncols_b])
    for i in range(0,numspecies_b):
        for m in range (0,numspecies_a):
            minputchar_b_mask[m,i]=lcompatibplantas[m]
    return minputchar_a_mask,minputchar_b_mask,lcompatibplantas    

def bino_mutual(filename,year_periods,hay_foodweb,hay_superpredadores,data_save='',dirtrabajo='',direntrada='',dirsal='',\
                eliminarenlaces=0,pl_ext=[],pol_ext=[],os='',fichreport='',com='', algorithm='MoMutualism', plants_blossom_prob=1.0,\
                plants_blossom_sd=0.01, plants_blossom_type = 'Binary', blossom_pert_list='', verbose=True,exit_on_extinction=False,\
                N0plants='',N0pols='',release='',Bssvar_period=0.1,Bssvar_sd=0.0,Bssvar_modulationtype_list=[],Bssvar_species=[]):
    global hayreport
    global frep
    global output_verbose
    global cuentaperpl,cuentaperpol
    global Logistic_abs
    global model_r_alpha
    global pendiente,blossomperiod,sd,periodo
    
    systemextinction = False
    periods = year_periods * diasdelanio
    May = (algorithm=='May')
    haymut = (algorithm!='NoMutualism')    
    model_r_alpha =  (algorithm=='Verhulst') or (algorithm=='NoMutualism')
    #print("model_ra"+str(model_r_alpha))
    Logistic_abs = (algorithm=='Logistic_abs')   
    hayreport=len(fichreport)>0
    output_verbose = verbose
    if hayreport: 
        frep=open(fichreport,'w',encoding='utf-8')
    informa ("Binomial simulated mutualistic interaction. Input file: %s" %(filename))
    informa ('============================================================================')
    if len(com)>0: informa("User Comment: %s" % com)
    tinic=time()
    days_year=diasdelanio
    informa('Span: %d years' % (year_periods))
    informa('ALGORITHM: '+algorithm)
    informa('Release '+ (("%.02f")%(release/100)))
    filename_a=filename+'_a.txt'
    dt = dirtrabajo.replace('\\','/')
    if hayreport: frep.write("<br>Plants matrix: <a href='file:///"+dt+"/input/"+filename_a+"' target=_BLANK>"+filename_a+"<a><br>")
    l_minputchar_a=dlmreadlike(filename_a,direntrada)
    ''' If N0plants provided by command line'''
    if len(N0plants)>0:
        l_minputchar_a[-5][0] = int(N0plants)
    minputchar_a = np.array(l_minputchar_a,dtype=float)
    try:
        nrows_a=len(minputchar_a)
        rdaux = float(minputchar_a[-1,0])
        if (rdaux < 0):
            raise
        if (rdaux < 0) or (rdaux > 5):
            raise
    except:
        print("INPUT FILE BAD FORMAT")
        return
    ncols_a=len(minputchar_a[0])
    numspecies_a = ncols_a
    rowNindividuals_a, Alpha_a, cAlpha_a, r_a, rd_a, Nindividuals_a, ra_eff, ra_equs = init_lists_pop(periods,numspecies_a)

    if (len(blossom_pert_list)>0):
        lcompatibplantas = calc_compatib_plantas(numspecies_a,plants_blossom_prob,blossom_pert_list)
    else:
        lcompatibplantas = calc_compatib_plantas(numspecies_a,plants_blossom_prob,[g for g in range(0,numspecies_a)])
    filename_b=filename+'_b.txt'
    if hayreport: frep.write("Pollinators matrix: <a href='file:///"+dt+"/input/"+filename_b+"' target=_BLANK>"+filename_b+"<a><br>")
    l_minputchar_b=dlmreadlike(filename_b,direntrada)
    ''' If N0pols provided by command line'''
    if len(N0pols)>0:
        l_minputchar_b[-5][0] = N0pols
    minputchar_b = np.array(l_minputchar_b,dtype=float)
    nrows_b=len(minputchar_b)
    ncols_b=len(minputchar_b[0])
    numspecies_b=ncols_b
    rowNindividuals_b, Alpha_b, cAlpha_b, r_b, rd_b, Nindividuals_b, rb_eff, rb_equs = init_lists_pop(periods,numspecies_b)

    K_c=[]
    Nindividuals_c=[]
    rowNindividuals_c=[]
    r_c=[]
    r_cperiod=[]
    minputchar_c=[]
    numspecies_c =  0
    
    if hay_foodweb>0:
        filename_c=filename+'_c.txt'
        filename_d=filename+'_d.txt'
        if hayreport: frep.write("Predators matrix c:<a href='file:///"+dt+"/input/"+filename_c+"' target=_BLANK>"+filename_c+"<a><br>")
        frep.write("Predators matrix d:<a href='file:///"+dt+"/input/"+filename_d+"' target=_BLANK>"+filename_d+"<a><br>")
        l_minputchar_c=dlmreadlike(filename_c,direntrada)
        minputchar_c = np.array(l_minputchar_c,dtype=float)
        nrows_c=len(minputchar_c)
        ncols_c=len(minputchar_c[0])
        numspecies_c=ncols_c;
        informa("Predator species : %d" %numspecies_c)
        for n in range(numspecies_c):
            rowNindividuals_c.append(int(minputchar_c[nrows_c-3][n]))
            K_c.append(int(minputchar_c[nrows_c-2][n]))
            r_c.append(minputchar_c[nrows_c-1][n])
            r_cperiod.append(calc_r_periodo_vh(minputchar_c[nrows_c-1][n],invperiod))
        Nindividuals_c.append(rowNindividuals_c)
        l_minputchar_d=dlmreadlike(filename_d,direntrada)
        minputchar_d = np.array(l_minputchar_d,dtype=float)
        nrows_d=len(minputchar_d)
        ncols_d=len(minputchar_d[0])
        numspecies_d = ncols_d
  
    informa("Plant species: %d" %numspecies_a)
    for n in range(numspecies_a):
        rowNindividuals_a.append(int(minputchar_a[-5][n]))
        cAlpha_a.append(float(minputchar_a[-4][n]))
        Alpha_a.append(float(minputchar_a[-3][n]))
        r_a.append(minputchar_a[-2][n])
        rd_a.append(minputchar_a[-1][n])
    Nindividuals_a[0]=np.array(rowNindividuals_a)
    
    informa("Plant initial populations %s" % rowNindividuals_a)
    if (plants_blossom_type=='Binary'):
        informa("Blossom probability %s, type %s. Plant affected species:%s" % (plants_blossom_prob,plants_blossom_type,str(blossom_pert_list)))
    else:
        informa("Blossom probability, type %s, mean %s, standard dev. %s. Plant affected species:%s" % (plants_blossom_type,plants_blossom_prob,plants_blossom_sd,str(blossom_pert_list)))
    informa("Pollinator species: %d" %numspecies_b)
    for n in range(numspecies_b):
        rowNindividuals_b.append(int(minputchar_b[-5][n]))
        cAlpha_b.append(float(minputchar_b[-4][n]))
        Alpha_b.append(float(minputchar_b[-3][n]))
        r_b.append(minputchar_b[-2][n])
        rd_b.append(minputchar_b[-1][n])
    Nindividuals_b[0]=np.array(rowNindividuals_b)
   
    informa("Pollinator initial populations %s" % rowNindividuals_b)
    if eliminarenlaces>0:
        cuenta=cuentaenlaces(minputchar_a)
        hayqueborrar=math.floor(eliminarenlaces*cuenta)
        informa ("Links %d. Will be deleted %d" % (cuentaenlaces(minputchar_a),hayqueborrar))
        if hayqueborrar > 0:
            periodoborr=periods/(hayqueborrar+1)
        else:
            periodoborr = periods

    # Extinction analysis
    inicioextplantas=inicioextpolin=cuentaperpl=cuentaperpol=0
    hayextplantas=len(pl_ext)>0
    hayextpolin = len(pol_ext) > 0
    j=0
    if hayextplantas and (pl_ext['species'][0]=='ALL'):
        pl_ext['species'] = []
        for m in range(0,numspecies_a+1):
            pl_ext['species'].append(m)
    if hayextpolin and (pol_ext['species'][0]=='ALL'):
        pol_ext['species'] = []
        for m in range(0,numspecies_a+1):
            pol_ext['species'].append(m)
    if (blossom_pert_list[0]=='ALL'):
        bloss_species = []
        for m in range(0,numspecies_a+1):
            bloss_species.append(m)
    else:
        bloss_species = blossom_pert_list[:]
    nperpl, inicioextplantas, periodoextpl, spikepl, nperpol,\
    inicioextpolin, periodoextpol, spikepol = init_perturbations(pl_ext, pol_ext, year_periods, inicioextplantas,\
                                                                 inicioextpolin, hayextplantas, hayextpolin)
    links_deletion = (eliminarenlaces>0)
    
    # Blossom variability analysis
    pBssvar_species = []
    hay_bssvar = (Bssvar_sd>0.0) & (len(Bssvar_species)>0)
    if (hay_bssvar):
        if (Bssvar_modulationtype_list[0]=='linear'):
            strvalor = (", Slope %0.2f " % (Bssvar_modulationtype_list[1]))
        else:
            if (Bssvar_modulationtype_list[0]=='sin'):
                strvalor = (", Period %0.2f " % (Bssvar_modulationtype_list[1]))
            else:
                strvalor=''
        straff = (". Affected species:%s" % (Bssvar_species))
        informa("Blossom variability active. Period: %0.2f (%% of year), Initial moment standard dev.: %0.04f (%% of year), Type: %s%s%s "\
                % (Bssvar_period,Bssvar_sd,Bssvar_modulationtype_list[0],strvalor,straff))
        if (str(Bssvar_species[0]).upper()=='ALL'):
            listaspecies = list(range(numspecies_a))
        else:
            listaspecies = Bssvar_species
        #print(listaspecies)
        bssvar_allones = np.ones((year_periods,),dtype=float)
        for i in range(numspecies_a):
            if i in listaspecies:
                if (Bssvar_modulationtype_list[0]=='None'):
                    indspecies, varspecies = calcbssvarespecie(Bssvar_period,Bssvar_sd,sinmodulacion,year_periods)
                elif (Bssvar_modulationtype_list[0]=='linear'):
                    pendiente = float(Bssvar_modulationtype_list[1])
                    indspecies, varspecies = calcbssvarespecie(Bssvar_period,Bssvar_sd,lineal,year_periods)
                elif (Bssvar_modulationtype_list[0]=='sin'):
                    periodo = int(Bssvar_modulationtype_list[1])
                    indspecies, varspecies = calcbssvarespecie(Bssvar_period,Bssvar_sd,sinusoidal,year_periods)            
                pBssvar_species.append(np.array(varspecies))
            else:
                pBssvar_species.append(bssvar_allones)
          
        #print(pBssvar_species)
         
    for k in range (periods-1):
        ''' The compatibilty matrixes masks are created when the year starts ''' 
        #if ((k%diasdelanio)==0):
        if not(k%diasdelanio):                      # Much faster than if ((k%diasdelanio)==0)    
            if (not(systemextinction)):
                if (algorithm!='Verhulst'):
                    if (k>diasdelanio) and ((np.array(lcompatibplantas)< plants_blossom_prob).sum()==0) and \
                       ((ra_equs[k-1]>0).sum()==0) and ((rb_equs[k-1]>0).sum()==0):
                        systemextinction = True
                        informa("ALARM !!!. System will collapse. Day %d (year %d)" % (k, k//diasdelanio))
                        if exit_on_extinction:
                            return(Nindividuals_a,Nindividuals_b,Nindividuals_c,ra_eff,rb_eff,ra_equs,ra_equs,0,0,0,0,0,0,systemextinction)
                else:
                    if (k>diasdelanio) and ((np.array(lcompatibplantas)< plants_blossom_prob).sum()==0) and \
                       ((ra_eff[k-1]>=tol_extincion).sum()==0) and ((rb_eff[k-1]>=tol_extincion).sum()==0):
                        systemextinction = True
                        informa("ALARM !!!. System will collapse. Day %d (year %d)" % (k, k//diasdelanio))
                        if exit_on_extinction:
                            return(Nindividuals_a,Nindividuals_b,Nindividuals_c,ra_eff,rb_eff,ra_equs,ra_equs,0,0,0,0,0,0,systemextinction)
            minputchar_a_mask, minputchar_b_mask, lcompatibplantas = calc_blossom_effect(numspecies_a,nrows_a,ncols_a,\
                                                                                         nrows_b,ncols_b,numspecies_b,plants_blossom_type,\
                                                                                         plants_blossom_prob,plants_blossom_sd,blossom_pert_list=bloss_species[:])
            minpeq_a = minputchar_a*minputchar_a_mask
            if hay_bssvar:
                for t in range(0,numspecies_b):
                    for l in range (0,numspecies_a):
                        minpeq_a[t,l] = minpeq_a[t,l]*pBssvar_species[l][k//diasdelanio]
            minpeq_b = minputchar_b*minputchar_b_mask
        # Eliminacion aleatoria de enlaces
        if (links_deletion):
            if (k>0) and (k%periodoborr==0):
                minputchar_a,fil,col=borraenlace(minputchar_a)
                minputchar_b[col][fil]=0
                minpeq_a = minputchar_a*minputchar_a_mask
                minpeq_b = minputchar_b*minputchar_b_mask
                informa("Day:%d (year %d). Deleted links plant %d <-> pollinator %d" % (k, k//diasdelanio ,fil,col))
        [populations_evolution(n,"Plant",numspecies_a,algorithm,hay_foodweb, pl_ext, May, haymut,\
                                             days_year, ra_eff, ra_equs, minpeq_a, j, cAlpha_a, Alpha_a, r_a, rd_a, Nindividuals_a,\
                                             numspecies_b, Nindividuals_b, Nindividuals_c, minputchar_c, numspecies_c,\
                                             inicioextplantas, hayextplantas, nperpl, periodoextpl, spikepl, k,\
                                             model_r_alpha) for n in range(numspecies_a)]
        ra_eff[0,]=ra_eff[1,]
        [populations_evolution(n,"Pollinator",numspecies_b,algorithm,hay_foodweb, pol_ext, May, haymut,\
                                              days_year, rb_eff, rb_equs, minpeq_b, j, cAlpha_b, Alpha_b, r_b, rd_b, Nindividuals_b,\
                                              numspecies_a, Nindividuals_a, Nindividuals_c, minputchar_c, numspecies_c,\
                                              inicioextpolin, hayextpolin, nperpol, periodoextpol, spikepol,k,\
                                              model_r_alpha) for n in range(numspecies_b)]
        rb_eff[0,]=rb_eff[1,]
        if hay_foodweb:
            rowNi=[]
            for n in range (numspecies_c):
                coef_bij_matrix=0
                c_devorados=0
                r_eqsum=0
                signo=signfunc(r_c[n])   
                if (Nindividuals_c[k][n]>0):
                    for j in range(numspecies_a):
                        coef_bij_matrix=Nindividuals_a[k][j]*minputchar_d[j][n]
                        r_eqsum+=coef_bij_matrix
                    for j in range(numspecies_b):
                        coef_bij_matrix=Nindividuals_b[k][j]*minputchar_d[j+numspecies_a][n]
                        r_eqsum+=coef_bij_matrix
                    rperiodequivalente = calc_r_periodo_vh(r_eqsum,invperiod)
                    term_c=round(Nindividuals_c[k][n]*(1-Nindividuals_c[k][n]/K_c[n]))
                    if term_c>1:
                        incPredatoria=np.random.binomial(term_c,-math.expm1(-rperiodequivalente))
                    else:
                        incPredatoria=0
                    incNmalth= np.random.binomial(Nindividuals_c[k][n],-math.expm1(-calc_r_periodo_vh(abs(r_c[n]),invperiod)))
                    signo=signfunc(r_c[n])
                    pop_c=Nindividuals_c[k][n]+incPredatoria+signo*incNmalth #+ incNlogistic-c_devorados
                    if (pop_c==0):
                        informa("Predator species %d extinction in day %d" % (n,k))
                else:
                    pop_c=0;
                rowNi.append(pop_c)
            Nindividuals_c.append(rowNi)
        
    maxa_individuos = np.max(Nindividuals_a)
    maxb_individuos = np.max(Nindividuals_b)
    max_reff = max(np.max([ra_eff]),np.max([rb_eff]))
    min_reff = min(np.min([ra_eff]),np.min([rb_eff]))
    max_requs = max(np.max([ra_equs]),np.max([rb_equs]))
    min_requs = min(np.min([ra_equs]),np.min([rb_equs]))

    tfin=time()
    informa ("Elapsed time %.02f s" % (tfin-tinic))
    speriodos = str(int(periods/diasdelanio))
    if (data_save==1):
        nsal = dlmwritelike(filename+'_'+algorithm+"_a_populations_",speriodos,Nindividuals_a,dirsal,os)
        rsal = dlmwritelike(filename+'_'+algorithm+"_a_rs_",speriodos,ra_eff,dirsal,os)
        requsal = dlmwritelike(filename+'_'+algorithm+"_a_requs_",speriodos,ra_equs,dirsal,os)
        if hayreport: 
            frep.write("<br>Plant populations data: <a href='"+nsal+"' target=_BLANK'>"+nsal+"<a>")
            frep.write("<br>Plant effective rates data: <a href='"+rsal+"' target=_BLANK'>"+rsal+"<a>")
            frep.write("<br>Plant equivalent rates data: <a href='"+requsal+"' target=_BLANK'>"+requsal+"<a>")
        nsal=dlmwritelike(filename+'_'+algorithm+"_b_populations_",speriodos,Nindividuals_b,dirsal,os)
        rsal = dlmwritelike(filename+'_'+algorithm+"_b_rs_",speriodos,rb_eff,dirsal,os)
        requsal = dlmwritelike(filename+'_'+algorithm+"_b_requs_",speriodos,rb_equs,dirsal,os)
        if hayreport: 
            frep.write("<br>Pollinators evolution data: <a href='"+nsal+"' target=_BLANK'>"+nsal+"<a>")
            frep.write("<br>Pollinators effective rates data: <a href='"+rsal+"' target=_BLANK'>"+rsal+"<a>")
            frep.write("<br>Pollinators equivalent rates data: <a href='"+requsal+"' target=_BLANK'>"+requsal+"<a>")
        if hay_foodweb>0:
            nsal=dlmwritelike(filename+'_'+algorithm+"_c",speriodos,Nindividuals_c,dirsal,os)
            if hayreport: frep.write("Predators evolution data: <a href='"+nsal+"' target=_BLANK'>"+nsal+"<a><br>")
    informa ('')
    informa ('Created %s' % datetime.datetime.now())
    if hayreport: frep.close()
    return(Nindividuals_a,Nindividuals_b,Nindividuals_c,ra_eff,rb_eff,ra_equs,rb_equs,maxa_individuos,maxb_individuos,\
           max_reff,min_reff,max_requs,min_requs,systemextinction,pBssvar_species)

def phase_map(na,nb,speciesa,speciesb,ciclosinic,ciclosfin):
    auxa=[]
    auxb=[]
    for k in range (ciclosinic,ciclosfin):
        auxa.append(na[k][speciesa])
        auxb.append(nb[k][speciesb])
    plt.hexbin(auxa,auxb,bins='log',cmap=cm.YlOrRd)
    plt.show()
    return


def setxtickssubplot(displayinic, periods, a):
    a.set_xlim([displayinic,periods])
    if ((displayinic + periods) //365 >= 10):
        ninter = 10
    else:
        ninter = ((displayinic + periods) //365)
    
    intervalo = (displayinic + periods) / ninter
    rangodias = np.arange(displayinic, periods + diasdelanio, intervalo)
    rangoanios = rangodias / diasdelanio
    a.xaxis.set_ticks(rangodias)
    xlabels = list(rangoanios)
    a.set_xticklabels(xlabels)
    a.grid(True)


def pintasubplot(na, min_individuos, max_individuos, displayinic, periods, factorescala, numspecies, titulo, ylabel):
    global ax
    plt.title(titulo)
    plt.ylabel(ylabel)
    for i in range(numspecies):
        graf = []
        x = []
        for k in range(displayinic, periods):
            graf.append(na[k][i])
            x.append(k)
        
        plt.plot(x, graf, color=cm.Set1(i / (numspecies)), lw=2)
        ax.plot(0,0, color=cm.Set1(i/(numspecies)),label='%i' % i)
    a = plt.gca()
    a.set_ylim([0-factorescala*abs(min_individuos), factorescala * max_individuos])
    setxtickssubplot(displayinic, periods, a)
    if numspecies < 11:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fancybox=True, shadow=True)
        leg = plt.gca().get_legend()
        ltext  = leg.get_texts()
        plt.setp(ltext, fontsize='small')
    

def mutual_render(na,nb,ra_eff,rb_eff,ra_equs,rb_equs,maxa_individuos,maxb_individuos,max_reff,min_reff,max_equs,min_equs,\
                  filename,displayinic,periods,dirsalida,algorithm='',fichreport='',os='',dirtrabajo='',Bssvar_coefs=[]):
    global ax
    # Si los valores de reff son muy pequenios, se cambia la escala de esas graficas
    matplotlib.rc('xtick', labelsize=8) 
    matplotlib.rc('ytick', labelsize=8) 
    years = periods/diasdelanio
    hayreport=len(fichreport)>0
    if hayreport: 
        frep=open(fichreport,'a',encoding='utf-8')
    factorescala=1.1
    numspecies_a=len(na[0])
    numspecies_b=len(nb[0])
    plt.figure('Mutualist network simulation. Input file: '+filename,dpi=resolucion,figsize=(ancho,alto))
    
    
    ax=plt.subplot(3, 2, 1)
    pintasubplot(na, 0, maxa_individuos, displayinic, periods, factorescala, numspecies_a, 'Plants', 'Individuals')
    
    ax=plt.subplot(3, 2, 3)
    pintasubplot(ra_eff, min_reff, max_reff, displayinic, periods, factorescala, numspecies_a, '','Efficient growth rate')
    plt.xlabel('Years')
    
    
    if (len(Bssvar_coefs)):
        ax=plt.subplot(3, 2, 5)
        plt.ylabel('Blossom variability coeffs.')
        plt.xlabel('Years')
        plt.grid(True)
        for i in range(numspecies_a):
            graf=[]
            x=[]
            for k in range (0,numanyos):
                graf.append(Bssvar_coefs[i][k])
                x.append(1+k)
            plt.plot(x,graf,color=cm.Set1(i/(numspecies_a)),lw=2)
            ax.plot(0,0, color=cm.Set1(i/(numspecies_a)),label='%i' % i)
        a = plt.gca()
        a.set_ylim([0,1.1])
        if numspecies_a < 11:
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
            ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fancybox=True, shadow=True)
            leg = plt.gca().get_legend()
            ltext  = leg.get_texts()
            plt.setp(ltext, fontsize='small')
    
    '''
    ax=plt.subplot(3, 2, 5)
    plt.ylabel('Equivalent growth rate')
    plt.xlabel('Days')
    plt.grid(True)
    for i in range(numspecies_a):
        graf=[]
        x=[]
        for k in range (displayinic,periods-1):
            graf.append(ra_equs[k][i])
            x.append(k)
        plt.plot(x,graf,color=cm.Set1(i/(numspecies_a)),lw=2)
        ax.plot(0,0, color=cm.Set1(i/(numspecies_a)),label='%i' % i)
    a = plt.gca()
    a.set_ylim([0-factorescala*abs(min_equs),factorescala*max_equs])
    if (numspecies_a<11):
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fancybox=True, shadow=True)
        leg = plt.gca().get_legend()
        ltext  = leg.get_texts()
        plt.setp(ltext, fontsize='small')
    '''
    
    ax=plt.subplot(3, 2, 2)
    pintasubplot(nb, 0, maxb_individuos, displayinic, periods, factorescala, numspecies_b, 'Polllinators', '')
    ax=plt.subplot(3, 2, 4)
    pintasubplot(rb_eff, min_reff, max_reff, displayinic, periods, factorescala, numspecies_b, '','')
    plt.xlabel('Years')
    
    '''
    ax=plt.subplot(3, 2, 6)
    #plt.ylabel('Equivalent growth rate')
    plt.xlabel('Days')
    plt.grid(True)
    for i in range(numspecies_b):
        graf=[]
        x=[]
        for k in range (displayinic,periods-1):
            graf.append(rb_equs[k][i])
            x.append(k)
        plt.plot(x,graf,color=cm.Paired(i/(numspecies_b)),lw=2)
        ax.plot(0,0, color=cm.Paired(i/(numspecies_b)),label='%i' % i)
    a = plt.gca()
    a.set_ylim([0-factorescala*abs(min_equs),factorescala*max_equs])
    if numspecies_b<11:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fancybox=True, shadow=True)
        leg = plt.gca().get_legend()
        ltext  = leg.get_texts()
        plt.setp(ltext, fontsize='small')
     '''
    
    dt = dirtrabajo.replace('\\','/');    
    nsal= 'output_pict_plantsandpols_'+filename+'_'+algorithm+'_'+os+'_'+str(years)+'.png'
    plt.savefig(str(dt+'/'+dirsalida.replace('\\','/')+nsal),bbox_inches=0)
    if hayreport:
        frep.write("<p align=left><br>Populations evolution picture<br>")
        frep.write("<IMG SRC=file:///%s ALIGN=LEFT  width=1200 BORDER=0>" % str(dt+'/'+dirsalida.replace('\\','/')+nsal)) 
        frep.write('</p>')
        frep.write('<br>')
        frep.close()
    #plt.show(2)
    plt.close()
    
def food_render(na,nb,nc,maxa_individuos,maxb_individuos,filename,displayinic,periods,dirsalida,algorithm='',fichreport='',os='',dirtrabajo=''):
    hayreport=len(fichreport)>0
    if hayreport: 
        frep=open(fichreport,'a',encoding='utf-8')
    matplotlib.rc('xtick', labelsize=8) 
    matplotlib.rc('ytick', labelsize=8) 
    factorescala=1.2
    numspecies_a=len(na[0])
    numspecies_b=len(nb[0])
    numspecies_c=len(nc[0])
    plt.figure('Mutualist network simulation. Input file: '+filename,dpi=resolucion,figsize=(ancho,alto))
    ax=plt.subplot(3, 1, 1)
    plt.title('Plants')
    plt.ylabel('Individuals')
    #plt.xlabel('Days')
    plt.grid(True)
    for i in range(numspecies_a):
        graf=[]
        x=[]
        for k in range (displayinic,periods):
            graf.append(na[k][i])
            x.append(k)
        plt.plot(x,graf,color=cm.Set1(i/(numspecies_a)),lw=2)
    a = plt.gca()
    a.set_ylim([0,factorescala*maxa_individuos])
    if numspecies_b<11:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    ax=plt.subplot(3, 1, 2)
    plt.title('Pollinators')
    plt.ylabel('Individuals')
    #plt.xlabel('Days')
    plt.grid(True)
    for i in range(numspecies_b):
        graf=[]
        x=[]
        for k in range (displayinic,periods):
            graf.append(nb[k][i])
            x.append(k)
        plt.plot(x,graf,color=cm.Paired(i/(numspecies_b)),lw=2)
    a = plt.gca()
    a.set_ylim([0,factorescala*maxb_individuos])
    if numspecies_b<11:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    ax=plt.subplot(3, 1, 3)
    plt.title('Predators')
    plt.ylabel('Individuals')
    plt.xlabel('Years')
    plt.grid(True)
    lmaxc = max(nc)
    maxc_individuos = max(lmaxc)
    for i in range(numspecies_c):
        graf=[]
        x=[]
        for k in range (displayinic,periods):
            graf.append(nc[k][i])
            x.append(k)
        plt.plot(x,graf,color=cm.Paired(i/(numspecies_c)),lw=2)
    a = plt.gca()
    a.set_ylim([0,factorescala*maxc_individuos])   
    if numspecies_c<11:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    nsal= 'output_pict_foodweb_'+filename+'_'+algorithm+'_'+os+'_'+str(years)+'.png'
    dt = dirtrabajo.replace('\\','/');
    plt.savefig(str(dt+'/'+dirsalida.replace('\\','/')+nsal),bbox_inches=0)
    plt.close()
    if hayreport:
        frep.write("<P align=left><br>Foodweb effect picture<br>")
        dt = dirtrabajo.replace('\\','/');
        frep.write("<IMG SRC=file:///%s ALIGN=LEFT  width=1200 BORDER=0>" % str(dt+'/'+dirsalida.replace('\\','/')+nsal)) 
        frep.write('<p>')
        frep.close()
    #plt.show(2)    
            