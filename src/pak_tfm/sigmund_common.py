'''
Created on 04/09/2014

@author: algarra
'''

import sys
import numpy as np
import datetime
import copy
import sigmund_GLOBALS as sgGL


class SimulationConditions():
    def __init__(self, filename ='', year_periods = '', hay_foodweb = False, 
                hay_superpredadores = False ,
                data_save='', dirtrabajo ='', direntrada='', dirsal='',
                eliminarenlaces=0, pl_ext=[], pol_ext=[], os='', fichreport='',
                com='', algorithm='MoMutualism', plants_blossom_prob=1.0,
                plants_blossom_sd=0.01, plants_blossom_type='Binary', 
                blossom_pert_list='', verbose=True, exit_on_extinction=False,
                N0plants='', N0pols='', release='', Bssvar_period=0.1, 
                Bssvar_sd=0.0, Bssvar_modulationtype_list=[], Bssvar_species=[]):
        self.filename = filename
        self.year_periods = year_periods 
        self.hay_foodweb = hay_foodweb
        self.hay_superpredadores = hay_superpredadores
        self.data_save = data_save
        self.dirtrabajo = dirtrabajo
        self.direntrada = direntrada
        self.dirsal = dirsal
        self.eliminarenlaces = eliminarenlaces
        self.pl_ext = copy.deepcopy(pl_ext)
        self.pol_ext = copy.deepcopy(pol_ext)
        self.os = os
        self.fichreport = fichreport
        self.com = com
        self.algorithm = algorithm
        self.plants_blossom_prob = plants_blossom_prob
        self.plants_blossom_sd = plants_blossom_sd
        self.plants_blossom_type = plants_blossom_type
        self.blossom_pert_list = blossom_pert_list
        self.verbose = verbose
        self.exit_on_extinction = exit_on_extinction
        self.N0plants = N0plants
        self.N0pols = N0pols
        self.release = release
        self.Bssvar_period = Bssvar_period
        self.Bssvar_sd = Bssvar_sd
        self.Bssvar_modulationtype_list = copy.deepcopy(Bssvar_modulationtype_list)
        self.Bssvar_species = copy.deepcopy(Bssvar_species)

class CanalInfo():
    """ CanalInfo is a wrapper of status information channels such as stdout and
    the report file """
    def __init__(self, device, newlinestr):
        self.device = device
        self.newline = newlinestr
    def write(self, texto):
        self.device.write(texto + self.newline)
    def close(self):
        self.device.close()
        
def inform_user(canalesinfo, texto):
    """ (list of CanalInfo, str) -> NoneType

    Execution status information for human user
    
    >>> ldevices = []
    >>> ldevices.append(CanalInfo(sys.stdout,'\\n'))
    >>> inform_user(ldevices, 'Info for you')
    Info for you
    """
    [canal.write(texto) for canal in canalesinfo]
    
def close_info_channels(canalesinfo):
    [canal.close() for canal in canalesinfo]
    
def open_info_channels(verbose, fichreport, mode):
    ldev = []
    lfich = []
    if verbose:    
        sout = CanalInfo(sys.stdout, '\n')
        ldev.append(sout)
    if len(fichreport) > 0:
        frep = CanalInfo(open(fichreport, mode, encoding='utf-8'), '<br>')
        ldev.append(frep)
        lfich.append(frep)
    return ldev, lfich

def dlmreadlike(inputfile, direntrada):
    try:
        data = open(direntrada + inputfile)
        mtx_input = []
        for each_line in data:
            try:
                mtx_input.append(each_line.strip().split('\t'))
            except ValueError:
                pass
        data.close()
        return(mtx_input)
    except IOError:
        print('The datafile %s is missing!' % inputfile)
        return(0)
 
def dlmwritelike(inputfile, nperiod, Nin, dirsalida, os):
    dsal = dirsalida.replace('\\', '/')
    nsal = 'output_data_' + inputfile + '_' + os + '_' + str(nperiod) + '.txt'
    print ("Output file %s" % dsal + nsal)
    salida = open(dsal + nsal, 'w', encoding='utf-8')
    for linea in Nin:
        for i in range(len(linea) - 1):
            y = "%.012f" % linea[i]
            salida.write(y + ',')
        salida.write(str(linea[-1]) + '\n');
    salida.close()
    return(nsal)

def find_max_values(Nindividuals_a, Nindividuals_b, ra_eff, rb_eff, ra_equs, 
                    rb_equs):   
    maxa_individuos = np.max(Nindividuals_a)
    maxb_individuos = np.max(Nindividuals_b)
    max_reff = max(np.max([ra_eff]), np.max([rb_eff]))
    min_reff = min(np.min([ra_eff]), np.min([rb_eff]))
    max_requs = max(np.max([ra_equs]), np.max([rb_equs]))
    min_requs = min(np.min([ra_equs]), np.min([rb_equs]))
    return maxa_individuos, maxb_individuos, max_reff, min_reff,\
           max_requs, min_requs

def start_report(ldev_inf, filename, com, year_periods, algorithm, release):
    inform_user(ldev_inf, \
      "Binomial simulated mutualistic interaction. Input file: %s" % (filename))
    inform_user(ldev_inf, 60 * '=')
    if len(com) > 0: inform_user(sgGL.ldev_inf, "User Comment: %s" % com)
    inform_user(ldev_inf, 'Span: %d years' % (year_periods))
    inform_user(ldev_inf, 'ALGORITHM: ' + algorithm)
    inform_user(ldev_inf, 'Release ' + (("%.02f") % (release / 100)))

def end_report(ldev_inf, lfich_inf, tfin, tinic, periods, data_save, 
               filename, algorithm, dirsal, os, Nindividuals_a, ra_eff, ra_equs,
               Nindividuals_b, rb_eff, rb_equs, Nindividuals_c, hay_food_web):    
    inform_user(ldev_inf, "Elapsed time %.02f s" % (tfin - tinic))
    speriodos = str(int(periods / sgGL.DAYS_YEAR))
    if (data_save == 1):
        nsal = dlmwritelike(filename + '_' + algorithm + "_a_populations_", 
                            speriodos, Nindividuals_a, dirsal, os)
        rsal = dlmwritelike(filename + '_' + algorithm + "_a_rs_", speriodos, 
                            ra_eff, dirsal, os)
        requsal = dlmwritelike(filename + '_' + algorithm + "_a_requs_", 
                               speriodos, ra_equs, dirsal, os)
        inform_user(lfich_inf, "Plant populations data: <a href='"\
                          + nsal + "' target=_BLANK'>" + nsal + "<a>")
        inform_user(lfich_inf, "Plant effective rates data: <a href='"\
                          + rsal + "' target=_BLANK'>" + rsal + "<a>")
        inform_user(lfich_inf, "Plant equivalent rates data: <a href='"\
                          + requsal + "' target=_BLANK'>" + requsal + "<a>")
        nsal = dlmwritelike(filename + '_' + algorithm + "_b_populations_", 
                            speriodos, Nindividuals_b, dirsal, os)
        rsal = dlmwritelike(filename + '_' + algorithm + "_b_rs_", 
                            speriodos, rb_eff, dirsal, os)
        requsal = dlmwritelike(filename + '_' + algorithm + "_b_requs_", 
                               speriodos, rb_equs, dirsal, os)
        inform_user(lfich_inf, "Pollinators evolution data: <a href='"\
                          + nsal + "' target=_BLANK'>" + nsal + "<a>")
        inform_user(lfich_inf,
                          "Pollinators effective rates data: <a href='" +\
                           rsal + "' target=_BLANK'>" + rsal + "<a>")
        inform_user(lfich_inf, 
                          "Pollinators equivalent rates data: <a href='" +\
                          requsal + "' target=_BLANK'>" + requsal + "<a>")
        if hay_food_web:
            nsal = dlmwritelike(filename + '_' + algorithm + "_c", speriodos, 
                                Nindividuals_c, dirsal, os)
            inform_user(lfich_inf, 
                              "Predators evolution data: <a href='" \
                              + nsal + "' target=_BLANK'>" + nsal + "<a><br>")
    inform_user(ldev_inf, '')
    inform_user(ldev_inf, 'Created %s' % datetime.datetime.now())
    close_info_channels(lfich_inf)

if __name__ == '__main__':
    import doctest
    doctest.testmod()