'''
Created on 04/09/2014

@author: algarra
'''
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.backends.backend_tkagg
import os
import sys
import numpy as np

global ancho
global alto
global resolucion
global ax

import sigmund_GLOBALS as sgGL
import sigmund_common as sgcom

ancho = 16
alto = 10
resolucion = 600

def setxtickssubplot(displayinic, periods, a, periodsinyears = False):
    if (not(periodsinyears)):
        a.set_xlim([displayinic, periods])
        if ((displayinic + periods) // sgGL.DAYS_YEAR >= 10):
            ninter = 10
        else:
            ninter = ((displayinic + periods) // sgGL.DAYS_YEAR) 
        intervalo = (displayinic + periods) / ninter
        rangodias = np.arange(displayinic, periods + sgGL.DAYS_YEAR, intervalo)
        rangoanios = rangodias / sgGL.DAYS_YEAR
        a.xaxis.set_ticks(rangodias)
        xlabels = list(rangoanios)
    else:
        a.set_xlim([displayinic, periods])
        if (periods >= 10):
            ninter = 10
        else:
            ninter = periods
        intervalo = (displayinic + periods ) / ninter
        rangoanios = np.arange(displayinic, periods + 1, intervalo)
        a.xaxis.set_ticks(rangoanios)
        xlabels = []
        for i in range(0,len(rangoanios)):
            xlabels.append("%0.1f" % (rangoanios[i]))
    a.set_xticklabels(xlabels)
    a.grid(True)


def display_legend():
    global ax
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fancybox=True, 
              shadow=True)
    handles, labels = ax.get_legend_handles_labels()
    labels = [chr(ord(i) + 1) for i in labels]
    ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1, 1), 
              fancybox=True, shadow=True)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize='small')

def pintasubplot(na, min_value, max_value, displayinic, periods, 
                 factorescala, numspecies, titulo, ylabel, 
                 periodsinyears = False):
    global ax
    plt.title(titulo)
    plt.ylabel(ylabel)
#     if periodsinyears:
#         periods += 1
    for i in range(numspecies):
        graf = []
        x = []
        for k in range(displayinic, periods):
            graf.append(na[k][i])
            x.append(k)
        if periodsinyears:
            graf.append(na[-1][i])
            x.append(1+k)
        plt.plot(x, graf, color=cm.Set1(i / (numspecies)), 
                 lw=calc_lw_width(numspecies))
        ax.plot(0, 0, color=cm.Set1(i / (numspecies)), label='%i' % i)
    a = plt.gca()
    a.set_ylim([-0.01 - factorescala * abs(min_value),\
                factorescala * max_value])
    setxtickssubplot(displayinic, periods, a, periodsinyears)
    if numspecies < 11:
        display_legend()
 
def mutual_render(simulation_params, sig_ret_val, displayinic, periods,  
                  verbose=True):
    global ax
    if (len(sig_ret_val.pBssvar_species)):
        nrows = 3
    else:
        nrows = 2
    matplotlib.rc('xtick', labelsize=8) 
    matplotlib.rc('ytick', labelsize=8) 
    years = periods // sgGL.DAYS_YEAR
    sgGL.ldev_inf, sgGL.lfich_inf = sgcom.open_info_channels(verbose, 
                                              simulation_params.fichreport, 'a')    
    factorescala = 1.1
    numspecies_a = len(sig_ret_val.Nindividuals_a[0])
    numspecies_b = len(sig_ret_val.Nindividuals_b[0])
    plt.figure('Mutualist network simulation. Input file: ' +\
               simulation_params.filename, dpi=resolucion, 
               figsize=(ancho, alto))
    ax = plt.subplot(nrows, 2, 1)
    pintasubplot(sig_ret_val.Nindividuals_a, 0, sig_ret_val.maxa_individuos, 
                 displayinic, periods, factorescala, 
                 numspecies_a, 'Plants', 'Individuals')
    ax = plt.subplot(nrows, 2, 3)
    pintasubplot(sig_ret_val.ra_eff, sig_ret_val.min_reff, sig_ret_val.max_reff,
                 displayinic, periods, factorescala, numspecies_a, '', 
                 'Efficient growth rate')
    plt.xlabel('Years')
    if (len(sig_ret_val.pBssvar_species)):
        for i in range(numspecies_a):
            graf = [sig_ret_val.pBssvar_species[i][0]]
            x = [0]
            for k in range (0, years):
                graf.append(sig_ret_val.pBssvar_species[i][k])
                x.append(1 + k)
        ax = plt.subplot(nrows, 2, 5)
        listacoefs = np.array(sig_ret_val.pBssvar_species)
        listacoefs = np.c_[listacoefs[:,0],listacoefs]
        listacoefs = list(listacoefs.transpose())
        pintasubplot(listacoefs, 0, 1, displayinic, years, factorescala,
                     numspecies_a, '', 'Blossom variability coeffs.', 
                     periodsinyears = True) 
        plt.xlabel('Years')    
    ax = plt.subplot(nrows, 2, 2)
    pintasubplot(sig_ret_val.Nindividuals_b, 0, sig_ret_val.maxb_individuos, 
                 displayinic, periods, factorescala,
                 numspecies_b, 'Polllinators', '')
    ax = plt.subplot(nrows, 2, 4)
    pintasubplot(sig_ret_val.rb_eff, sig_ret_val.min_reff, sig_ret_val.max_reff,
                 displayinic, periods, factorescala,
                 numspecies_b, '', '')
    plt.xlabel('Years')
    dt = simulation_params.dirtrabajo.replace('\\', '/');    
    nsal = 'output_pict_plantsandpols_' + simulation_params.filename +\
           '_' + simulation_params.algorithm + '_' + simulation_params.os +\
           '_' + str(years) + '.png'
    plt.savefig(str(dt + '/' + simulation_params.dirsal.replace('\\', '/') + nsal), 
                bbox_inches=0)
    sgcom.inform_user(sgGL.lfich_inf, "<p align=left>Populations evolution picture")
    sgcom.inform_user(sgGL.lfich_inf, \
                      "<IMG SRC=file:///%s ALIGN=LEFT  width=1200 BORDER=0>" %\
                      str(dt + '/' + simulation_params.dirsal.replace('\\', '/') + nsal)) 
    sgcom.inform_user(sgGL.lfich_inf, '</p><br>')
    sgcom.close_info_channels(sgGL.lfich_inf)
    plt.close()
    
def calc_lw_width(numspecies):
    return(0.5)

def food_render(simulation_params, sig_ret_val, displayinic, periods, 
                verbose=True):
    global ax
    sgGL.ldev_inf, sgGL.lfich_inf = sgcom.open_info_channels(verbose, 
                                            simulation_params.fichreport, 'a')    
    
    matplotlib.rc('xtick', labelsize=8) 
    matplotlib.rc('ytick', labelsize=8) 
    factorescala = 1.2
    numspecies_a = len(sig_ret_val.Nindividuals_a[0])
    numspecies_b = len(sig_ret_val.Nindividuals_b[0])
    numspecies_c = len(sig_ret_val.Nindividuals_c[0])
    plt.figure('Mutualist network simulation. Input file: ' +\
               simulation_params.filename,
               dpi=resolucion, figsize=(ancho, alto))
    ax = plt.subplot(3, 1, 1)
    plt.title('Plants')
    plt.ylabel('Individuals')
    # plt.xlabel('Days')
    plt.grid(True)
    for i in range(numspecies_a):
        graf = []
        x = []
        for k in range (displayinic, periods):
            graf.append(sig_ret_val.Nindividuals_a[k][i])
            x.append(k)
        plt.plot(x, graf, color=cm.Set1(i / (numspecies_a)),
                 lw=calc_lw_width(numspecies_a))
    a = plt.gca()
    a.set_ylim([0, factorescala * sig_ret_val.maxa_individuos])
    if numspecies_b < 11:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    ax = plt.subplot(3, 1, 2)
    plt.title('Pollinators')
    plt.ylabel('Individuals')
    # plt.xlabel('Days')
    plt.grid(True)
    for i in range(numspecies_b):
        graf = []
        x = []
        for k in range (displayinic, periods):
            graf.append(sig_ret_val.Nindividuals_b[k][i])
            x.append(k)
        plt.plot(x, graf, color=cm.Paired(i / (numspecies_b)),
                 lw=calc_lw_width(numspecies_b))
    a = plt.gca()
    a.set_ylim([0, factorescala * sig_ret_val.maxb_individuos])
    if numspecies_b < 11:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    ax = plt.subplot(3, 1, 3)
    plt.title('Predators')
    plt.ylabel('Individuals')
    plt.xlabel('Years')
    plt.grid(True)
    lmaxc = max(sig_ret_val.Nindividuals_c)
    maxc_individuos = max(lmaxc)
    for i in range(numspecies_c):
        graf = []
        x = []
        for k in range (displayinic, periods):
            graf.append(sig_ret_val.Nindividuals_c[k][i])
            x.append(k)
        plt.plot(x, graf, color=cm.Paired(i / (numspecies_c)),
                 lw=calc_lw_width(numspecies_c))
    a = plt.gca()
    a.set_ylim([0, factorescala * maxc_individuos])   
    if numspecies_c < 11:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    nsal = 'output_pict_foodweb_' + simulation_params.filename + '_' +\
            simulation_params.algorithm + '_' + simulation_params.os +\
           '_' + str(periods) + '.png'
    dt = simulation_params.dirtrabajo.replace('\\', '/');
    plt.savefig(str(dt + '/' + simulation_params.dirsal.replace('\\', '/') + nsal), 
                bbox_inches=0)
    plt.close()
    sgcom.inform_user(sgGL.lfich_inf, "<P align=left><br>Foodweb effect picture<br>")
#    dt = simulation_params.dirtrabajo.replace('\\', '/');
    sgcom.inform_user(sgGL.lfich_inf,\
                      "<IMG SRC=file:///%s ALIGN=LEFT  width=1200 BORDER=0>" %\
                      str(dt + '/' + simulation_params.dirsal.replace('\\', '/') + nsal)) 
    sgcom.inform_user(sgGL.lfich_inf, '<p>')
    sgcom.close_info_channels(sgGL.lfich_inf)   
     
if __name__ == '__main__':
    pass