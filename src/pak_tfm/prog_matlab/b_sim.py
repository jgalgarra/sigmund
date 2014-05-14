'''
Created on 10/05/2012

@author: Javier Garcia Algarra

This module includes the mutualist algorithm
'''
from math import exp
from time import time
import datetime
import numpy as np
import sys
from numpy.random import binomial
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.backends.backend_tkagg
import os

global ancho
global alto
global resolucion
global pathout

ancho=13
alto=8
resolucion=600

def cuentaenlaces(mat_in):
    cuenta=0;
    fil=len(mat_in)-3
    col=len(mat_in[0])
    for i in range(fil):
        for j in range(col):
            if abs(mat_in[i][j])>0:
                cuenta+=1
    return cuenta

def borraenlace(mat_in):
    intentos=0
    fil=len(mat_in)-3
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

def dlmwritelike(inputfile,nperiod,Nin,dirsalida,os):
    nsal= 'output_data_'+inputfile+'_'+os+'_'+str(nperiod)+'.txt'
    print ("Output file %s" %dirsalida+nsal)
    salida=open(dirsalida+nsal,'w',encoding='utf-8')
    for linea in Nin:
        for i in range(len(linea)-1):
            salida.write(str(linea[i])+'\t')
        salida.write(str(linea[-1])+'\n');
    salida.close()
    return(nsal)


def r_periodo(Ranual,period):
    #if Ranual>=1:
    #    return(Ranual/period)
    return ((1+Ranual)**(1/period)-1)
    #return (Ranual/period)

def informa(text):
    print(text,file=sys.stdout)
    if hayreport: frep.write(text+'<br>')

def CalcCoefMay(minputchar,r,K):
    return(minputchar if r!=0 else 0)
    #return(abs(K*minputchar/r) if r!=0 else 0)

def CicloMay(r_species, period, inctermMay, Nindivs, K):
    rspneq = r_periodo(abs(r_species),period)
    signosp = np.sign(r_species)
    termEq = Nindivs * (signosp-(Nindivs/K))
    incEq = binomial(round(abs(termEq)),1-exp(-1*rspneq)) if (abs(termEq)> 1) else 0
    #return(Nindivs -inc1Eq + inc2Eq)
    return(Nindivs + np.sign(termEq)*incEq + inctermMay)
    #return(Nindivs + np.sign(term3Eq)*inc3Eq)

def ValMutMay(r_species,beta,period,N1,N2,K1):
    rspneq = r_periodo(abs(r_species),period)
    termEq = round(beta*N1*N2/K1)
    incEq = binomial(round(abs(termEq)),1-exp(-1*rspneq)) if (abs(termEq)> 1) else 0
    return(incEq)

def CicloNewModel(rtot_species, period, Nindivs, K):
    rspneq = r_periodo(abs(rtot_species),period)
    term = round(np.sign(rtot_species)*Nindivs-Nindivs**2/K)
    if (term == 0):
      incNmalth = 0
    else:
      incNmalth= binomial(abs(term),1-exp(-1*rspneq)) 
    inc_pop= Nindivs+ np.sign(term)*incNmalth
    return(inc_pop)
                   
    
def bino_mutual(filename,periods,hay_foodweb,hay_superpredadores,data_save='',direntrada='',dirsal='',eliminarenlaces=0,pl_ext=[],pol_ext=[],os='',fichreport='',com='', algorithm='MoMutualism'):          
    global hayreport
    global frep
    May = (algorithm=='May')
    haymut = (algorithm!='NoMutualism')  
            
    hayreport=len(fichreport)>0
    if hayreport: 
        frep=open(fichreport,'w',encoding='utf-8')
    informa ("Binomial simulated mutualistic interaction. Input file: %s   days: %ld" %(filename,periods))
    informa ('============================================================================')
    if len(com)>0: informa("User Comment: %s" % com)
    tinic=time()
    period_year=365
    informa('ALGORITHM: '+algorithm)
    filename_a=filename+'_a.txt'
    print("filename_a "+filename_a)
    if hayreport: frep.write("<br>Plants matrix: <a href='../input/"+filename_a+"' target=_BLANK>"+filename_a+"<a><br>")
    minputchar_a=dlmreadlike(filename_a,direntrada)
    try:
        nrows_a=len(minputchar_a)
    except:
        print("INPUT FILE BAD FORMAT")
        return
    ncols_a=len(minputchar_a[0])
    for i in range (nrows_a):
        for j in range (ncols_a):
            minputchar_a[i][j]=float(minputchar_a[i][j])
    numspecies_a=ncols_a
    K_a=[]
    Nindividuals_a=[]
    rowNindividuals_a=[]
    r_a=[]
    rd_a=[]
    ra_eq=[]
    rb_eq=[]
    maxima_ind=[]
    minima_ind=[]
    maximb_ind=[]
    minimb_ind=[]
    maxim_req=[]
    minim_req=[]
    filename_b=filename+'_b.txt'
    if hayreport: frep.write("Pollinators matrix: <a href='../input/"+filename_b+"' target=_BLANK>"+filename_b+"<a><br>")
    minputchar_b=dlmreadlike(filename_b,direntrada)
    nrows_b=len(minputchar_b)
    ncols_b=len(minputchar_b[0])
    for i in range (nrows_b):
        for j in range (ncols_b):
            minputchar_b[i][j]=float(minputchar_b[i][j])
    numspecies_b=ncols_b;
    K_b=[]
    Nindividuals_b=[]
    rowNindividuals_b=[]
    r_b=[]
    rd_b=[]
    
    K_c=[]
    Nindividuals_c=[]
    rowNindividuals_c=[]
    r_c=[]
    r_cperiod=[]
    
    if hay_foodweb>0:
        filename_c=filename+'_c.txt'
        filename_d=filename+'_d.txt'
        if hayreport: frep.write("Predators matrix c:<a href='../input/"+filename_c+" target=_BLANK'>"+filename_c+"<a><br>")
        frep.write("Predators matrix d:<a href='../input/"+filename_d+" target=_BLANK'>"+filename_d+"<a><br>")
        minputchar_c=dlmreadlike(filename_c,direntrada)
        nrows_c=len(minputchar_c)
        ncols_c=len(minputchar_c[0])
        for i in range (nrows_c):
            for j in range (ncols_c):
                minputchar_c[i][j]=float(minputchar_c[i][j])
        numspecies_c=ncols_c;
        informa("Predator species : %d" %numspecies_c)
        for n in range(numspecies_c):
            rowNindividuals_c.append(int(minputchar_c[nrows_c-3][n]))
            K_c.append(int(minputchar_c[nrows_c-2][n]))
            r_c.append(minputchar_c[nrows_c-1][n])
            r_cperiod.append(r_periodo(minputchar_c[nrows_c-1][n],period_year))
        Nindividuals_c.append(rowNindividuals_c)
        minputchar_d=dlmreadlike(filename_d,direntrada)
        nrows_d=len(minputchar_d)
        ncols_d=len(minputchar_d[0])
        for i in range (nrows_d):
            for j in range (ncols_d):
                minputchar_d[i][j]=float(minputchar_d[i][j])
  
    informa("Plant species: %d" %numspecies_a)
    for n in range(numspecies_a):
        rowNindividuals_a.append(int(minputchar_a[-4][n]))
        K_a.append(int(minputchar_a[-3][n]))
        r_a.append(minputchar_a[-2][n])
        rd_a.append(minputchar_a[-1][n])
    Nindividuals_a.append(rowNindividuals_a)
    informa("Plant initial populations %s" % rowNindividuals_a)
    informa("Pollinator species: %d" %numspecies_b)
    for n in range(numspecies_b):
        rowNindividuals_b.append(int(minputchar_b[-4][n]))
        K_b.append(int(minputchar_b[-3][n]))
        r_b.append(minputchar_b[-2][n])
        rd_b.append(minputchar_b[-1][n])
    Nindividuals_b.append(rowNindividuals_b)
    informa("Pollinator initial populations %s" % rowNindividuals_b)
    if eliminarenlaces>0:
        cuenta=cuentaenlaces(minputchar_a)
        hayqueborrar=round(eliminarenlaces*cuenta)
        informa ("Links %d. Will be deleted %d" % (cuentaenlaces(minputchar_a),hayqueborrar))
        periodoborr=round(0.75*(periods-1)/(hayqueborrar+1));

    # Analisis de extinciones
    inicioextplantas=inicioextpolin=cuentaperpl=cuentaperpol=0
    hayextplantas=len(pl_ext)>0
    if hayextplantas:
        inicioextplantas=round(pl_ext['start']*periods)
        nperpl=pl_ext['numperiod']
        periodoextpl=pl_ext['period']
        spikepl=round(periodoextpl*pl_ext['spike'])
        print(spikepl)
        informa("Forced extinctions. Plants species %s, period (years): %d, numperiods: %d, spike: %0.2f, rate: %.03f, start: %.d" % (pl_ext['species'],periodoextpl/period_year,nperpl,pl_ext['spike'],float(pl_ext['rate']),inicioextplantas))
    hayextpolin=len(pol_ext)>0
    if hayextpolin:
        inicioextpolin=round(pol_ext['start']*periods)
        nperpol=pol_ext['numperiod']
        periodoextpol=pol_ext['period']
        spikepol=round(periodoextpol*pol_ext['spike'])
        informa("Forced extinctions. Pollinators species %s, period (years): %d, numperiods: %d, spike: %0.2f, rate: %.03f, start: %.d" % (pol_ext['species'],periodoextpol/period_year,nperpol,pol_ext['spike'],float(pol_ext['rate']),inicioextpolin))
    for k in range (periods-1):
        rowNi=[]
        row_req=[]
        # Eliminacion aleatoria de enlaces
        if (k>0) and (eliminarenlaces>0) and (k%periodoborr==0):
            minputchar_a,fil,col=borraenlace(minputchar_a)
            minputchar_b[col][fil]=0
            informa("Day %d. Deleted link %d,%d" % (k,fil,col))
        for n in range (numspecies_a):
            coef_bij_matrix=0;  
            r_eqsum=  0 
            term_May = 0
            a_devorados=0
            rtot_a = 0
            signo=np.sign(r_a[n])
            r_muerte= rd_a[n]
            if (Nindividuals_a[k][n]>0):
                # Extinciones de plantas
                nohuboextincion=1
                if hayextplantas and (cuentaperpl<nperpl):
                    if (k>=inicioextplantas) and (n in pl_ext['species']):
                        posindextpl=(k-inicioextplantas)%periodoextpl
                        if (posindextpl>=0) and (posindextpl< spikepl):
                            #r_eqsum=pl_ext['rate']
                            r_muerte = pl_ext['rate']
                            signo=-1
                            nohuboextincion=0
                        if (posindextpl==spikepl-1):
                            cuentaperpl+=1
                if nohuboextincion:
                    for j in range(numspecies_b):
                        if haymut:
                            #coef_bij_matrix=minputchar_a[j][n]*Nindividuals_b[k][j]/K_a[n]
                            coef_bij_matrix=minputchar_a[j][n]*Nindividuals_b[k][j]
                            r_eqsum+=coef_bij_matrix
                            term_May+=ValMutMay(r_a[n],CalcCoefMay(minputchar_a[j][n],r_a[n],K_a[n]),period_year,Nindividuals_a[k][n],Nindividuals_b[k][j],K_a[n])
                rtot_a = r_a[n]+r_eqsum-r_muerte
                # Efecto de los depredadores
                if hay_foodweb>0:
                    for j in range(numspecies_c):
                        if hay_foodweb>0:
                            if (minputchar_c[n][j]>0):
                                rceq=r_periodo(Nindividuals_c[k][j]*minputchar_c[n][j],period_year)
                                a_devorados+=binomial(Nindividuals_a[k][n],1-exp(-1*rceq))
                            #if (a_devorados>0):
                            #    print (k,n,a_devorados)               
                #incNmalth=signo*binomial(Nindividuals_a[k][n],1-exp(-1*rperiodequivalente))
                # New algorithm
                if not(May):
                    pop_a = CicloNewModel(rtot_a, period_year, Nindividuals_a[k][n] ,K_a[n])-a_devorados
                # May's classical model
                else: 
                    pop_a = CicloMay(r_a[n], period_year, term_May, Nindividuals_a[k][n], K_a[n])-a_devorados
                if (pop_a==0):
                    informa("Plant species %d extinction in day %d" % (n,k))
            else:
                pop_a=0
            rowNi.append(pop_a)
            #print(n,r_eqsum)
            row_req.append(rtot_a-(abs(rtot_a)*Nindividuals_a[k][n]/K_a[n]))
            #suma_a[n]+=pop_a
        Nindividuals_a.append(rowNi)
        ra_eq.append(row_req)
        maxima_ind.append(max(rowNi))
        minima_ind.append(min(rowNi))
        maxim_req.append(max(row_req))
        minim_req.append(min(row_req))
           
        rowNi=[]
        row_req=[]
        
        for n in range (numspecies_b):
            coef_bij_matrix=0
            b_devorados=0
            r_eqsum= 0 
            rtot_b = 0
            term_May = 0
            signo=np.sign(r_b[n])
            r_muerte = rd_b[n]   
            if (Nindividuals_b[k][n]>0):
                # Extinciones de bichos
                nohuboextincion=1
                if hayextpolin and (cuentaperpol<nperpol):
                    if (k>=inicioextpolin) and (n in pol_ext['species']):
                        posindextpol=(k-inicioextpolin)%periodoextpol
                        if (posindextpol>=0) and (posindextpol< spikepol):
                        #if ((k-inicioextpolin)%periodoextpol>=0) and ((k-inicioextpolin)%periodoextpol<spikepol):
                                #r_eqsum=pol_ext['rate']
                                r_muerte = pol_ext['rate']
                                signo=-1
                                nohuboextincion=0
                        if (posindextpol==spikepol-1):
                            cuentaperpol+=1
                if nohuboextincion:
                    for j in range(numspecies_a):
                        if haymut:
                            #coef_bij_matrix=minputchar_b[j][n]*Nindividuals_a[k][j]/K_b[n]
                            coef_bij_matrix=minputchar_b[j][n]*Nindividuals_a[k][j]
                            r_eqsum+=coef_bij_matrix
                            term_May+=ValMutMay(r_b[n],CalcCoefMay(minputchar_b[j][n],r_b[n],K_b[n]),period_year,Nindividuals_b[k][n],Nindividuals_a[k][j],K_b[n])
                rtot_b = r_b[n]+r_eqsum-r_muerte
                # Efecto de los depredadores
                if hay_foodweb>0:
                    for j in range(numspecies_c):
                        if hay_foodweb>0:
                            if (minputchar_c[n][j]>0):
                                rceq=r_periodo(Nindividuals_c[k][j]*minputchar_c[n][j],period_year)
                                b_devorados+=binomial(Nindividuals_b[k][n],1-exp(-1*rceq))
                if not(May):
                    pop_b = CicloNewModel(rtot_b, period_year, Nindividuals_b[k][n] ,K_b[n])-b_devorados
                # May's classical model
                else:
                    pop_b = CicloMay(r_b[n], period_year, term_May, Nindividuals_b[k][n], K_b[n])-b_devorados
                if (pop_b==0):
                    informa("Pollinator species %d extinction in day %d" % (n,k))
            else:
                pop_b=0
            rowNi.append(pop_b)
            row_req.append(rtot_b-(abs(rtot_b)*Nindividuals_b[k][n]/K_b[n]))
        Nindividuals_b.append(rowNi)
        rb_eq.append(row_req)
        maximb_ind.append(max(rowNi))
        minimb_ind.append(min(rowNi))
        maxim_req.append(max(row_req))
        minim_req.append(min(row_req))
        
        rowNi=[]
        if hay_foodweb>0:
            for n in range (numspecies_c):
                coef_bij_matrix=0
                c_devorados=0
                r_eqsum=0
                signo=np.sign(r_c[n])   
                if (Nindividuals_c[k][n]>0):
                    for j in range(numspecies_a):
                        #coef_bij_matrix=minputchar_b[j][n]*Nindividuals_a[k][j]/K_b[n]
                        coef_bij_matrix=Nindividuals_a[k][j]*minputchar_d[j][n]
                        r_eqsum+=coef_bij_matrix
                    for j in range(numspecies_b):
                        #coef_bij_matrix=minputchar_b[j][n]*Nindividuals_a[k][j]/K_b[n]
                        coef_bij_matrix=Nindividuals_b[k][j]*minputchar_d[j+numspecies_a][n]
                        r_eqsum+=coef_bij_matrix
                    rperiodequivalente=r_periodo(r_eqsum,period_year)
                    term_c=round(Nindividuals_c[k][n]*(1-Nindividuals_c[k][n]/K_c[n]))
                    #term_c=round(Nindividuals_c[k][n])
                    if term_c>1:
                        incPredatoria=binomial(term_c,1-exp(-1*rperiodequivalente))
                    else:
                        incPredatoria=0
                    incNmalth= binomial(Nindividuals_c[k][n],1-exp(-1*r_periodo(abs(r_c[n]),period_year)))
                    signo=np.sign(r_c[n])
                    pop_c=Nindividuals_c[k][n]+incPredatoria+signo*incNmalth #+ incNlogistic-c_devorados
                    if (pop_c==0):
                        informa("Predator species %d extinction in day %d" % (n,k))
                else:
                    pop_c=0;
                rowNi.append(pop_c)
            Nindividuals_c.append(rowNi)
        
    maxa_individuos=(max(maxima_ind))
    maxb_individuos=(max(maximb_ind))
    max_req=(max(maxim_req))
    min_req=(min(minim_req))
    tfin=time()
    informa ("Elapsed time %.02f s" % (tfin-tinic))
    if (data_save==1):
        nsal=dlmwritelike(filename+algorithm+"_a",periods,Nindividuals_a,dirsal,os)
        if hayreport: frep.write("<br>Plants evolution data: <a href='"+nsal+"' target=_BLANK'>"+nsal+"<a>")
        nsal=dlmwritelike(filename+algorithm+"_b",periods,Nindividuals_b,dirsal,os)
        if hayreport: frep.write("<br>Pollinators evolution data: <a href='"+nsal+"' target=_BLANK'>"+nsal+"<a><br>")
        if hay_foodweb>0:
            nsal=dlmwritelike(filename+algorithm+"_c",periods,Nindividuals_c,dirsal,os)
            if hayreport: frep.write("Predators evolution data: <a href='"+nsal+"' target=_BLANK'>"+nsal+"<a><br>")
    informa ('')
    informa ('Created %s' % datetime.datetime.now())
    if hayreport: frep.close()
    return(Nindividuals_a,Nindividuals_b,Nindividuals_c,ra_eq,rb_eq,maxa_individuos,maxb_individuos,max_req,min_req)

def phase_map(na,nb,speciesa,speciesb,ciclosinic,ciclosfin):
    auxa=[]
    auxb=[]
    for k in range (ciclosinic,ciclosfin):
        auxa.append(na[k][speciesa])
        auxb.append(nb[k][speciesb])
    plt.hexbin(auxa,auxb,bins='log',cmap=cm.YlOrRd)
    #plt.scatter(auxa,auxb,s=1,lw=1,marker='.',facecolor='0.5')
    plt.show()
    return

def mutual_render(na,nb,ra_eq,rb_eq,maxa_individuos,maxb_individuos,max_req,min_req,filename,displayinic,periods,dirsalida,fichreport='',os=''):
    # Si los valores de req son muy pequenios, se cambia la escala de esas graficas
    if (max_req-min_req<0.0001):
        max_req=0.01
        min_req=-0.01
    matplotlib.rc('xtick', labelsize=8) 
    matplotlib.rc('ytick', labelsize=8) 
    hayreport=len(fichreport)>0
    if hayreport: 
        frep=open(fichreport,'a',encoding='utf-8')
    factorescala=1.1
    numspecies_a=len(na[0])
    numspecies_b=len(nb[0])
    plt.figure('Mutualist network simulation. Input file: '+filename,dpi=resolucion,figsize=(ancho,alto))
    ax=plt.subplot(2, 2, 1)
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
        plt.plot(x,graf,color=plt.cm.Set1(i/(numspecies_a)),lw=1)
    a = plt.gca()
    a.set_ylim([0,factorescala*maxa_individuos])
    if numspecies_b<11:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    ax=plt.subplot(2, 2, 3)
    plt.ylabel('Equivalent growth rate')
    plt.xlabel('Days')
    plt.grid(True)
    for i in range(numspecies_a):
        graf=[]
        x=[]
        for k in range (displayinic,periods-1):
            graf.append(ra_eq[k][i])
            x.append(k)
        plt.plot(x,graf,color=plt.cm.Set1(i/(numspecies_a)),lw=1)
        ax.plot(0,0, color=plt.cm.Set1(i/(numspecies_a)),label='%i' % i)
    a = plt.gca()
    a.set_ylim([0-factorescala*abs(min_req),factorescala*max_req])
    if (numspecies_a<11):
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fancybox=True, shadow=True)
        leg = plt.gca().get_legend()
        ltext  = leg.get_texts()
        plt.setp(ltext, fontsize='small')
    #plt.plot(ra_eq,linewidth=0.2)
    #plt.show()
    ax=plt.subplot(2, 2, 2)
    plt.title('Pollinators')
    #plt.ylabel('Individuals')
    #plt.xlabel('Days')
    plt.grid(True)
    for i in range(numspecies_b):
        graf=[]
        x=[]
        for k in range (displayinic,periods):
            graf.append(nb[k][i])
            x.append(k)
        plt.plot(x,graf,color=plt.cm.Paired(i/(numspecies_b)),lw=1)
    a = plt.gca()
    a.set_ylim([0,factorescala*maxb_individuos])
    if numspecies_b<11:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    ax=plt.subplot(2, 2, 4)
    #plt.ylabel('Equivalent growth rate')
    plt.xlabel('Days')
    plt.grid(True)
    for i in range(numspecies_b):
        graf=[]
        x=[]
        for k in range (displayinic,periods-1):
            graf.append(rb_eq[k][i])
            x.append(k)
        plt.plot(x,graf,color=plt.cm.Paired(i/(numspecies_b)),lw=1)
        ax.plot(0,0, color=plt.cm.Paired(i/(numspecies_b)),label='%i' % i)
    a = plt.gca()
    a.set_ylim([0-factorescala*abs(min_req),factorescala*max_req])
    #plt.plot(rb_eq,linewidth=0.2)
    if numspecies_b<11:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fancybox=True, shadow=True)
        leg = plt.gca().get_legend()
        ltext  = leg.get_texts()
        plt.setp(ltext, fontsize='small')
    nsal= 'output_pict_plantsandpols_'+filename+'_'+os+'_'+str(periods)+'.png'
    plt.savefig((dirsalida+nsal),bbox_inches=0)
    plt.close()
    if hayreport:
        frep.write("<br><br>Populations evolution picture<br><table border=0><tr><td>")
        frep.write("<IMG SRC='%s' ALIGN=LEFT  width=900 BORDER=0>" % nsal) 
        frep.write('</td></tr></table>')
        frep.write('<P>')
        frep.close()
    #plt.show(2)
    
def food_render(na,nb,nc,maxa_individuos,maxb_individuos,filename,displayinic,periods,dirsalida,fichreport='',os=''):
    hayreport=len(fichreport)>0
    if hayreport: 
        frep=open(fichreport,'a',encoding='utf-8')
    matplotlib.rc('xtick', labelsize=8) 
    matplotlib.rc('ytick', labelsize=8) 
    factorescala=1.1
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
        plt.plot(x,graf,color=plt.cm.Set1(i/(numspecies_a)),lw=1)
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
        plt.plot(x,graf,color=plt.cm.Paired(i/(numspecies_b)),lw=1)
    a = plt.gca()
    a.set_ylim([0,factorescala*maxb_individuos])
    if numspecies_b<11:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    ax=plt.subplot(3, 1, 3)
    plt.title('Predators')
    plt.ylabel('Individuals')
    plt.xlabel('Days')
    plt.grid(True)
    for i in range(numspecies_c):
        graf=[]
        x=[]
        for k in range (displayinic,periods):
            graf.append(nc[k][i])
            x.append(k)
        plt.plot(x,graf,color=plt.cm.Paired(i/(numspecies_c)),lw=1)
    a = plt.gca()
    if (maxb_individuos>maxa_individuos):
        auxmaxc=maxb_individuos
    else:
        auxmaxc=maxa_individuos
    a.set_ylim([0,factorescala*auxmaxc])
    if numspecies_c<11:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    nsal= 'output_pict_foodweb_'+filename+'_'+os+'_'+str(periods)+'.png'
    plt.savefig((dirsalida+nsal),bbox_inches=0)
    plt.close()
    if hayreport:
        frep.write("<br><br>Foodweb effect picture<br>")
        frep.write("<IMG SRC='%s' NAME='graphics1' width=900 BORDER=0>" % nsal) 
        frep.write('<br>')
        frep.close()
    #plt.show(2)    
    
    
    
def dif_simulation(filename,periods):
    print ('Diferential simulation %s %ld' %(filename,periods))
    tinic=time()
    minputchar=dlmreadlike(filename)
    nrows=len(minputchar)
    ncols=len(minputchar[0])
    for i in range (nrows):
        for j in range (ncols):
            minputchar[i][j]=float(minputchar[i][j])
    numspecies=ncols
    print ("numspecies %d" %numspecies)
    K=[]
    Nindividuals=[]
    rowNindividuals=[]
    for n in range(numspecies):
        rowNindividuals.append(int(minputchar[numspecies][n]))
        K.append(int(minputchar[numspecies+1][n]))
    Nindividuals.append(rowNindividuals)
    period_year=365
    for k in range (periods-1):
        rowNi=[]
        for n in range (numspecies):
            r=minputchar[n][n]
            rperiod=r/period_year
            incNmalth= rperiod*Nindividuals[k][n]
            incNlogistic= rperiod*(Nindividuals[k][n]**2)/K[n]
            incNOtherspecies=0;
            rowNi.append(round(Nindividuals[k][n]+incNmalth-incNlogistic+incNOtherspecies))
        Nindividuals.append(rowNi)
    tfin=time()
    print ("Elapsed time %f s" % (tfin-tinic))
    dlmwritelike(filename,periods,Nindividuals,os)
    pobj=plt.plot(Nindividuals)
    return(pobj)



# bino_simulation('inputdata_2_muchos.txt',20000)