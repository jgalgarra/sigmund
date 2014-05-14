'''
Created on 10/05/2012

@author: remoto
'''
from math import *
from time import *
from numpy.random import binomial
import matplotlib.pyplot as plt


def dlmreadlike(inputfile):
    try:
        data = open(inputfile)
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

def dlmwritelike(inputfile,nperiod,Nin,method):
    nsal= 'output_periods_'+str(nperiod)+'_'+method+'_'+inputfile;
    print ("Output file %s" %nsal)
    salida=open(nsal,'w',encoding='utf-8')
    for linea in Nin:
        for i in range(len(linea)-1):
            salida.write(str(linea[i])+'\t')
        salida.write(str(linea[-1])+'\n');
    salida.close()


def bino_simulation(filename,periods):
    print ('Binomial simulation %s %ld' %(filename,periods))
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
            # Variation due to malthusian parameter r
            incNmalth=binomial(Nindividuals[k][n],1-exp(-1*rperiod));
            # Second term of logistic equation
            incNlogistic= binomial((Nindividuals[k][n]**2)/K[n],1-exp(-1*rperiod));
            # Terms due to other species
            incNOtherspecies=0;
            for j in range(numspecies):
               if (j!=n):
                  incNOtherspecies=incNOtherspecies+binomial(round(Nindividuals[k][j]*Nindividuals[k][n]/K[n]),1-exp(-1*rperiod*minputchar[n][j]));
            rowNi.append(round(Nindividuals[k][n]+incNmalth-incNlogistic+incNOtherspecies))  
        Nindividuals.append(rowNi)
    tfin=time()
    print ("Elapsed time %f s" % (tfin-tinic))
    dlmwritelike(filename,periods,Nindividuals,'bino')
    plt.plot(Nindividuals)
    plt.show()

def bino_mutual(filename,periods):
    print ('Binomial simulation of mutualistic interaction %s %ld' %(filename,periods))
    tinic=time()
    
    filename_a=filename+'_a.txt'
    minputchar_a=dlmreadlike(filename_a)
    nrows_a=len(minputchar_a)
    ncols_a=len(minputchar_a[0])
    for i in range (nrows_a):
        for j in range (ncols_a):
            minputchar_a[i][j]=float(minputchar_a[i][j])
    numspecies_a=ncols_a
    print ("numspecies a %d" %numspecies_a)
    K_a=[]
    Nindividuals_a=[]
    rowNindividuals_a=[]
    r_a=[]
    
    filename_b=filename+'_b.txt'
    minputchar_b=dlmreadlike(filename_b)
    nrows_b=len(minputchar_b)
    ncols_b=len(minputchar_b[0])
    for i in range (nrows_b):
        for j in range (ncols_b):
            minputchar_a[i][j]=float(minputchar_b[i][j])
    numspecies_b=nrows_b-3;
    print ("numspecies b %d" %numspecies_b)
    K_b=[]
    Nindividuals_b=[]
    rowNindividuals_b=[]
    r_b=[]
    
    for n in range(numspecies_a):
        rowNindividuals_a.append(int(minputchar_a[nrows_a-3][n]))
        K_a.append(int(minputchar_a[nrows_a-2][n]))
        r_a.append(minputchar_a[nrows_a-1][n])
    Nindividuals_a.append(rowNindividuals_a)
    period_year=365
    for k in range (periods-1):
        rowNi=[]
        for n in range (numspecies_a):
            rperiod=float(r_a[n]/period_year)
            # Variation due to malthusian parameter r
            incNmalth=binomial(Nindividuals_a[k][n],1-exp(-1*rperiod));
            # Second term of logistic equation
            incNlogistic= binomial((Nindividuals_a[k][n]**2)/K_a[n],1-exp(-1*rperiod));
            # Terms due to other species
            incNOtherspecies=0;
            for j in range(numspecies_a):
                incNOtherspecies=incNOtherspecies+binomial(round(Nindividuals_a[k][j]*Nindividuals_a[k][n]/K_a[n]),1-exp(-1*rperiod*minputchar_a[n][j]));
            rowNi.append(round(Nindividuals_a[k][n]+incNmalth-incNlogistic+incNOtherspecies))  
        Nindividuals_a.append(rowNi)
    tfin=time()
    print ("Elapsed time %f s" % (tfin-tinic))
    dlmwritelike(filename_a,periods,Nindividuals_a,'bino')
    plt.plot(Nindividuals_a)
    plt.show()
    
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
    dlmwritelike(filename,periods,Nindividuals,'dif')
    pobj=plt.plot(Nindividuals)
    return(pobj)



# bino_simulation('inputdata_2_muchos.txt',20000)