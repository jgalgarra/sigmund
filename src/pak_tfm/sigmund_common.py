'''
Created on 04/09/2014

@author: algarra
'''

import sys
import numpy as np
import datetime
import copy
import pickle
import sigmund_GLOBALS as sgGL

class SimulationConditions():
    def __init__(self, filename ='', year_periods = '', hay_foodweb = False, 
                hay_superpredadores = False ,
                data_save='', dirtrabajo ='', direntrada='', dirsal='',
                eliminarenlaces=0, pl_ext=[], pol_ext=[], output_suff ='', 
                fichreport = '', com = '', algorithm = 'MoMutualism', 
                plants_blossom_prob = 1.0, plants_blossom_sd = 0.01, 
                plants_blossom_type = 'Binary', blossom_pert_list = '', 
                verbose = True, exit_on_extinction = False,
                N0plants='', N0pols='', release='', Bssvar_data = ''):
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
        self.output_suff = output_suff
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
        self.Bssvar_data = Bssvar_data
    
    def write2file(self,filesim):
        with open(filesim, 'wb') as outfile:
            pickle.dump(self, outfile, pickle.HIGHEST_PROTOCOL)
         
class MaximaValues():
    def __init__(self,maxa_individuos, maxb_individuos, max_reff, min_reff,\
    max_requs, min_requs):
        self.maxa_individuos = maxa_individuos
        self.maxb_individuos = maxb_individuos
        self.max_reff = max_reff
        self.min_reff = min_reff
        self.max_requs = max_requs
        self.min_requs = min_requs
        
class BlossomVariability():
    def __init__(self, Bssvar_period, Bssvar_sd, Bssvar_modulationtype_list,
                 Bssvar_species):
        self.Bssvar_period = Bssvar_period
        self.Bssvar_sd = Bssvar_sd
        self.Bssvar_modulationtype_list = Bssvar_modulationtype_list
        self.Bssvar_species = Bssvar_species
        
class SimulationReturnValues():
    def __init__(self, Nindividuals_a, Nindividuals_b, Nindividuals_c, ra_eff, 
                 rb_eff, ra_equs, rb_equs, mval, systemextinction, 
                 pBssvar_species):
        self.systemextinction = systemextinction
        self.Nindividuals_a = Nindividuals_a
        self.Nindividuals_b = Nindividuals_b
        self.Nindividuals_c = Nindividuals_c
        self.ra_eff = copy.deepcopy(ra_eff)
        self.rb_eff = copy.deepcopy(rb_eff)
        self.ra_equs = copy.deepcopy(ra_equs)
        self.rb_equs = copy.deepcopy(ra_equs)
        self.maxminval = mval
        self.pBssvar_species = copy.deepcopy(pBssvar_species)

class ExternalPerturbationConditions():
    def __init__(self,nperpl, inicioextplantas, periodoextpl, spikepl, 
                 nperpol, inicioextpol, periodoextpol, spikepol):
        self.nperpl = nperpl
        self.inicioextplantas = inicioextplantas
        self.periodoextpl = periodoextpl
        self.spikepl = spikepl
        self.nperpol = nperpol
        self.inicioextpol = inicioextpol
        self.periodoextpol = periodoextpol
        self.spikepol = spikepol
        
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
        return(-1)
 
def dlmwritelike(input_file,sim_cond, nperiod, Nin):
    dsal = sim_cond.dirsal.replace('\\', '/')
    nsal = 'output_data_' + input_file + '_' + sim_cond.output_suff + '_'+\
           str(nperiod) + '.txt'
    print ("Output file %s" % dsal + nsal)
    salida = open(dsal + nsal, 'w', encoding='utf-8')
    for linea in Nin:
        for i in range(len(linea) - 1):
            y = "%.012f" % linea[i]
            salida.write(y + ',')
        salida.write(str(linea[-1]) + '\n');
    salida.close()
    return(nsal)

def read_simulation_matrix(filename, dirtrabajo, direntrada, str_guild, 
                           name_guild, N0_guild, lfich_inf = ''):    
    filename_x = filename + str_guild
    dt = dirtrabajo.replace('\\', '/')
    inform_user(lfich_inf, name_guild + " matrix: <a href='file:///" + dt +\
                      "/input/" + filename_x + "' target=_BLANK>" +\
                      filename_x + "<a>")
    l_minputchar_x = dlmreadlike(filename_x, direntrada)
    ''' If N0_guild provided by command line'''
    if len(N0_guild) > 0:
        l_minputchar_x[-sgGL.LINES_INFO_MATRIX][0] = int(N0_guild)
    minputchar_x = np.array(l_minputchar_x, dtype=float)
    try:
        nrows_x = len(minputchar_x)
    except:
        print("INPUT FILE BAD FORMAT")
    ncols_x = len(minputchar_x[0])
    numspecies_x = ncols_x
    return numspecies_x, minputchar_x, nrows_x, ncols_x

def find_max_values(Nindividuals_a, Nindividuals_b, ra_eff, rb_eff, ra_equs, 
                    rb_equs):   
    mval = MaximaValues(np.max(Nindividuals_a),np.max(Nindividuals_b),
                         max(np.max([ra_eff]), np.max([rb_eff])),
                         min(np.min([ra_eff]), np.min([rb_eff])),
                         max(np.max([ra_equs]), np.max([rb_equs])),
                         min(np.min([ra_equs]), np.min([rb_equs])) )
    return(mval)


def create_list_species_affected(speciestext):
    if (str(speciestext) == 'ALL'):
        return(list(['ALL']))
    auxspec = speciestext.split(',')
    auxlspec = []
    for numspe in auxspec:
        if len(numspe) == 1:
            auxlspec.append([int(numspe)])
        else:
            rango = numspe.split(':')
            auxlspec.append([j for j in range(int(rango[0]), 1+int(rango[-1]))])
    listb = []
    for subla in auxlspec:
        for j in list(subla):
            listb.append(j)    
    listb = sorted(set(listb))   
    return listb

def start_report(ldev_inf, filename, com, year_periods, algorithm, release, 
                 hay_foodweb):
    inform_user(ldev_inf, \
      "Binomial simulated mutualistic interaction. Input file: %s" % (filename))
    inform_user(ldev_inf, 60 * '=')
    if len(com) > 0: inform_user(sgGL.ldev_inf, "User Comment: %s" % com)
    inform_user(ldev_inf, 'Span: %d years' % (year_periods))
    inform_user(ldev_inf, 'ALGORITHM: ' + algorithm)
    inform_user(ldev_inf, 'Release ' + (("%.02f") % (release / 100)))
    if hay_foodweb:
        inform_user(ldev_inf, 'Food web superimposed') 

def create_file_suffix(algorithm,output_suffix,ciclos):
    return('_'.join([algorithm,output_suffix,str(int(ciclos))]) )

def create_fichreport_name(reportpath,input_fname,file_suffix):
    return reportpath + 'rep_' + input_fname +'_'+file_suffix+ '.html'

def create_results_filename(sim_cond,string_file):
    return sim_cond.filename + '_' + sim_cond.algorithm + string_file

def end_report(ldev_inf, lfich_inf, sim_cond, tfin, tinic, periods, 
               Nindividuals_a, ra_eff, ra_equs,
               Nindividuals_b, rb_eff, rb_equs, Nindividuals_c):    
    inform_user(ldev_inf, "Elapsed time %.02f s" % (tfin - tinic))
    speriodos = str(int(periods / sgGL.DAYS_YEAR))
    if (sim_cond.data_save == 1):
        nsal = dlmwritelike(create_results_filename(sim_cond,"_a_populations"), 
                            sim_cond,speriodos, Nindividuals_a)
        rsal = dlmwritelike(create_results_filename(sim_cond,"_a_rs"),
                            sim_cond,speriodos, ra_eff)
        requsal = dlmwritelike(create_results_filename(sim_cond,"_a_requs"), 
                               sim_cond,speriodos, ra_equs)
        inform_user(lfich_inf, "Plant populations data: <a href='"\
                          + nsal + "' target=_BLANK'>" + nsal + "<a>")
        inform_user(lfich_inf, "Plant effective rates data: <a href='"\
                          + rsal + "' target=_BLANK'>" + rsal + "<a>")
        inform_user(lfich_inf, "Plant equivalent rates data: <a href='"\
                          + requsal + "' target=_BLANK'>" + requsal + "<a>")
        nsal = dlmwritelike(create_results_filename(sim_cond,"_b_populations"), 
                            sim_cond,speriodos, Nindividuals_b)
        rsal = dlmwritelike(create_results_filename(sim_cond,"_b_rs"), 
                            sim_cond,speriodos, rb_eff)
        requsal = dlmwritelike(create_results_filename(sim_cond,"_b_requs"), 
                               sim_cond,speriodos, rb_equs)
        inform_user(lfich_inf, "Pollinators evolution data: <a href='"\
                          + nsal + "' target=_BLANK'>" + nsal + "<a>")
        inform_user(lfich_inf, "Pollinators effective rates data: <a href='" +\
                           rsal + "' target=_BLANK'>" + rsal + "<a>")
        inform_user(lfich_inf, "Pollinators equivalent rates data: <a href='" +\
                          requsal + "' target=_BLANK'>" + requsal + "<a>")
        if sim_cond.hay_foodweb:
            nsal = dlmwritelike(create_results_filename(sim_cond,"_c"), 
                                sim_cond,speriodos, Nindividuals_c)
            inform_user(lfich_inf, 
                              "Predators evolution data: <a href='" \
                              + nsal + "' target=_BLANK'>" + nsal + "<a><br>")
    inform_user(ldev_inf, '')
    inform_user(ldev_inf, 'Created %s' % datetime.datetime.now())
    close_info_channels(lfich_inf)

if __name__ == '__main__':
    import doctest
    doctest.testmod()