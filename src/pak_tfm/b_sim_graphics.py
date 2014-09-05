'''
Created on 03/09/2014

@author: algarra
'''

import numpy as np
import b_sim_GLOBALS as bsG
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib

global ax
global ancho
global alto
global resolucion

ancho = 16
alto = 10
resolucion = 600

     
def calc_lw_width(numspecies):
    return(0.5)

def setxtickssubplot(displayinic, periods, a):
    a.set_xlim([displayinic, periods])
    if ((displayinic + periods) // 365 >= 10):
        ninter = 10
    else:
        ninter = ((displayinic + periods) // 365)
    
    intervalo = (displayinic + periods) / ninter
    rangodias = np.arange(displayinic, periods + bsG.DAYS_YEAR, intervalo)
    rangoanios = rangodias / bsG.DAYS_YEAR
    a.xaxis.set_ticks(rangodias)
    xlabels = list(rangoanios)
    a.set_xticklabels(xlabels)
    a.grid(True)
    
def pintasubplot(na, min_individuos, max_individuos, displayinic, periods, 
                 factorescala, numspecies, titulo, ylabel):
    global ax
    plt.title(titulo)
    plt.ylabel(ylabel)
    for i in range(numspecies):
        graf = []
        x = []
        for k in range(displayinic, periods):
            graf.append(na[k][i])
            x.append(k)
         
        plt.plot(x, graf, color=cm.Set1(i / (numspecies)), 
                 lw=calc_lw_width(numspecies))
        ax.plot(0, 0, color=cm.Set1(i / (numspecies)), label='%i' % i)
    a = plt.gca()
    a.set_ylim([0 - factorescala * abs(min_individuos),\
                factorescala * max_individuos])
    setxtickssubplot(displayinic, periods, a)
    if numspecies < 11:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1), 
                  fancybox=True, shadow=True)
        leg = plt.gca().get_legend()
        ltext = leg.get_texts()
        plt.setp(ltext, fontsize='small')
        
def mutual_render(na, nb, ra_eff, rb_eff, ra_equs, rb_equs, maxa_individuos,
                  maxb_individuos, max_reff, min_reff, max_equs, min_equs, 
                  filename, displayinic, periods, dirsalida,  
                  algorithm='',
                  fichreport='', verbose=True, os='', dirtrabajo='',
                  Bssvar_coefs=[]):
    global ax
    # Si los valores de reff son muy pequenios, se cambia la escala de esas graficas
    matplotlib.rc('xtick', labelsize=8) 
    matplotlib.rc('ytick', labelsize=8) 
    years = periods / bsG.DAYS_YEAR
    #ldevices_info, lfich_info = open_info_channels(verbose, fichreport, 'a')        
    factorescala = 1.1
    numspecies_a = len(na[0])
    numspecies_b = len(nb[0])
    plt.figure('Mutualist network simulation. Input file: ' +\
               filename, dpi=resolucion, figsize=(ancho, alto))
    ax = plt.subplot(3, 2, 1)
    pintasubplot(na, 0, maxa_individuos, displayinic, periods, factorescala, 
                 numspecies_a, 'Plants', 'Individuals')
    ax = plt.subplot(3, 2, 3)
    pintasubplot(ra_eff, min_reff, max_reff, displayinic, periods, factorescala,
                 numspecies_a, '', 'Efficient growth rate')
    plt.xlabel('Years')
     
    if (len(Bssvar_coefs)):
        ax = plt.subplot(3, 2, 5)
        plt.ylabel('Blossom variability coeffs.')
        plt.xlabel('Years')
        plt.grid(True)
        for i in range(numspecies_a):
            graf = []
            x = []
            for k in range (0, numanyos):
                graf.append(Bssvar_coefs[i][k])
                x.append(1 + k)
            plt.plot(x, graf, color=cm.Set1(i / (numspecies_a)),
                     lw=calc_lw_width(numspecies_a))
            ax.plot(0, 0, color=cm.Set1(i / (numspecies_a)), label='%i' % i)
        a = plt.gca()
        a.set_ylim([0, 1.1])
        if numspecies_a < 11:
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
            ax.legend(loc='upper left', bbox_to_anchor=(1, 1),
                      fancybox=True, shadow=True)
            leg = plt.gca().get_legend()
            ltext = leg.get_texts()
            plt.setp(ltext, fontsize='small')
         
    ax = plt.subplot(3, 2, 2)
    pintasubplot(nb, 0, maxb_individuos, displayinic, periods, factorescala,
                 numspecies_b, 'Polllinators', '')
    ax = plt.subplot(3, 2, 4)
    pintasubplot(rb_eff, min_reff, max_reff, displayinic, periods, factorescala,
                 numspecies_b, '', '')
    plt.xlabel('Years')
       
    dt = dirtrabajo.replace('\\', '/');    
    nsal = 'output_pict_plantsandpols_' + filename + '_' + algorithm + '_' +\
            os + '_' + str(years) + '.png'
    plt.savefig(str(dt + '/' + dirsalida.replace('\\', '/') + nsal), bbox_inches=0)
#     show_info_to_user(lfich_info, "<p align=left>Populations evolution picture")
#     show_info_to_user(lfich_info, \
#                       "<IMG SRC=file:///%s ALIGN=LEFT  width=1200 BORDER=0>" %\
#                       str(dt + '/' + dirsalida.replace('\\', '/') + nsal)) 
#     show_info_to_user(lfich_info, '</p><br>')
#     close_info_channels(lfich_info)
    plt.close()
     