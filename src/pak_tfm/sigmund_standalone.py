
import argparse
import pickle
import b_sim
import os
import sys
import sigmund_GLOBALS as sgGL
import sigmund_release as sgRL
import sigmund_common as sgcom
import sigmund_graphs as sggraph
import sgtools_common as sgtcom


def modify_float_par(condition,threshold,value):
    if float(condition)>threshold:
        return float(condition)
    else:
        return value
    
parser = argparse.ArgumentParser()

parser.add_argument('-simfile', action='store', dest='simfile_name',
                    help='Simulation parameters file')

parser.add_argument('-g', action='store_true', default=False,
                    dest='gengraphs',
                    help='Generate and store graphs')

parser.add_argument('-v', action='store_true', default=False,
                    dest='verbose',
                    help='Verbose output')

parser.add_argument('-years', action='store', dest='sim_years', default = '',
                    help='Simulation span in years')

parser.add_argument('-fw', action='store_true', default=False,
                    dest='foodweb',
                    help='Superimposed food web')

parser.add_argument('-outsf', action='store', dest='outsf', default='',
                    help='Append suffix to output file names')

parser.add_argument('-ds', action='store_true', default=False,
                    help='Results data save')

parser.add_argument('-dsdir', action='store', dest='dsdir', default='',
                    help='Data save directory')

parser.add_argument('-stop', action='store_true', default=False,
                    dest='stop_extinction',
                    help='On extinction stop simulation')
 
parser.add_argument('-rlink', action='store', dest='rlink', default = -1.0,
                    help='Randomlinks removal')

parser.add_argument('-Blprob', action='store', dest='Blprob', default = 0.0,
                     help='Blossom probability')

parser.add_argument('-Blsd', action='store', dest='Blsd', default = 0.0,
                     help='Blossom deviation')

parser.add_argument('-Bltype', action='store', dest='Bltype', default = '',
                     help='Blossom type: Binary, Gaussian')

parser.add_argument('-Blspecies',action='store', dest='Blspecies', default = '',
                     help='Blossom species')

parser.add_argument('-Bssvarper', action='store', dest='Bssvarper', default = 0.0,
                     help='Blossom variability period')

parser.add_argument('-Bssvarsd', action='store', dest='Bssvarsd', default = 0.0,
                     help='Blossom variability deviation')

parser.add_argument('-Bssvartype', action='store', dest='Bssvartype', default = '',
                     help='Blossom variability modulation type: None, linear, sin')

parser.add_argument('-Bssvartype_param', action='store', dest='Bssvartype_param', default = 0.0,
                     help='Blossom variability modulation parameter')

parser.add_argument('-Bssvarspecies',action='store', dest='Bssvarspecies', default = '',
                     help='Blossom variability affected species')

parser.add_argument('-pl_ext_period', action='store', dest='pl_ext_period', default = '',
                     help='Plants external perturbation, period in years')

parser.add_argument('-pl_ext_spike', action='store', dest='pl_ext_spike', default = '',
                     help='Plants external perturbation, spike (fraction of period)')

parser.add_argument('-pl_ext_numperiod', action='store', dest='pl_ext_numperiod', default = '',
                     help='Plants external perturbation, repeat perturbation')

parser.add_argument('-pl_ext_rate', action='store', dest='pl_ext_rate', default = '',
                     help='Plants external perturbation, rate')

parser.add_argument('-pl_ext_start', action='store', dest='pl_ext_start', default = '',
                     help='Plants external perturbation, initial year')

parser.add_argument('-pl_ext_species', action='store', dest='pl_ext_species', default = '',
                     help='Plants external perturbation, affected species')

parser.add_argument('-pol_ext_period', action='store', dest='pol_ext_period', default = '',
                     help='Pollinators external perturbation, period in years')

parser.add_argument('-pol_ext_spike', action='store', dest='pol_ext_spike', default = '',
                     help='Pollinators external perturbation, spike (fraction of period)')

parser.add_argument('-pol_ext_numperiod', action='store', dest='pol_ext_numperiod', default = '',
                     help='Pollinators external perturbation, repeat perturbation')

parser.add_argument('-pol_ext_rate', action='store', dest='pol_ext_rate', default = '',
                     help='Pollinators external perturbation, rate')

parser.add_argument('-pol_ext_start', action='store', dest='pol_ext_start', default = '',
                     help='Pollinators external perturbation, initial year')

parser.add_argument('-pol_ext_species', action='store', dest='pol_ext_species', default = '',
                     help='Pollinators external perturbation, affected species')

    
conditions = parser.parse_args()

spars = sgtcom.read_file_input(conditions.simfile_name)

filename = spars.filename
year_periods = spars.year_periods 
hay_foodweb = spars.hay_foodweb 
hay_superpredadores = spars.hay_superpredadores
data_save = spars.data_save
dirtrabajo = spars.dirtrabajo 
direntrada = spars.direntrada
dirsal = spars.dirsal
eliminarenlaces = spars.eliminarenlaces
pl_ext = spars.pl_ext
pol_ext = spars.pol_ext 
output_suff = spars.output_suff
fichreport = spars.fichreport.replace('\\', '/')
com = spars.com 
algorithm = spars.algorithm
plants_blossom_prob = spars.plants_blossom_prob 
plants_blossom_sd = spars.plants_blossom_sd
plants_blossom_type = spars.plants_blossom_type 
blossom_pert_list = spars.blossom_pert_list
release=spars.release
Bssvar_data = spars.Bssvar_data

if len(conditions.sim_years)>0:
    year_periods = int(conditions.sim_years)
if len(conditions.dsdir)>0:
    dirsal = dirsal + conditions.dsdir +'/'
    try:
        os.stat(dirsal)
    except:
        os.makedirs(dirsal)
verbose = conditions.verbose
data_save = conditions.ds
hay_foodweb = conditions.foodweb
eliminarenlaces = modify_float_par(conditions.rlink, -1.0, eliminarenlaces) 
plants_blossom_prob = modify_float_par(conditions.Blprob,0.0,plants_blossom_prob)
plants_blossom_sd = modify_float_par(conditions.Blsd,0.0,plants_blossom_sd)
if conditions.Bltype != '':
    plants_blossom_type = conditions.Bltype
if len(conditions.Blspecies)>0:
    blossom_pert_list = sgcom.create_list_species_affected(conditions.Blspecies)
Bssvar_data.Bssvar_period = modify_float_par(conditions.Bssvarper,0.0,Bssvar_data.Bssvar_period)
Bssvar_data.Bssvar_sd = modify_float_par(conditions.Bssvarsd,0.0,Bssvar_data.Bssvar_sd)
if conditions.Bssvartype != '':
    Bssvar_data.Bssvar_modulationtype_list[0] = conditions.Bssvartype
if float(conditions.Bssvartype_param)>0.0:
    if len(Bssvar_data.Bssvar_modulationtype_list)>1:
        Bssvar_data.Bssvar_modulationtype_list[1] = float(conditions.Bssvartype_param)
    else:
        Bssvar_data.Bssvar_modulationtype_list.append(float(conditions.Bssvartype_param))
if len(conditions.Bssvarspecies)>0:
    Bssvar_data.Bssvar_species = sgcom.create_list_species_affected(conditions.Bssvarspecies)
if len(conditions.pl_ext_period)>0:
    pl_ext['period'] = int(conditions.pl_ext_period)*sgGL.DAYS_YEAR
if len(conditions.pl_ext_numperiod)>0:
    pl_ext['numperiod'] = int(conditions.pl_ext_numperiod)
if len(conditions.pl_ext_spike)>0:
    pl_ext['spike'] = float(conditions.pl_ext_spike)
if len(conditions.pl_ext_rate)>0:
    pl_ext['rate'] = float(conditions.pl_ext_rate)
if len(conditions.pl_ext_start)>0:
    pl_ext['start'] = int(conditions.pl_ext_start)
if len(conditions.pl_ext_species)>0:
    if (conditions.pl_ext_species=='ALL'):
        pl_ext['species'] = ['ALL']
    else:
        pl_ext['species'] = [int(i) for i in conditions.pl_ext_species.split(',')]

if len(conditions.pol_ext_period)>0:
    pol_ext['period'] = int(conditions.pol_ext_period)*sgGL.DAYS_YEAR
if len(conditions.pol_ext_numperiod)>0:
    pol_ext['numperiod'] = int(conditions.pol_ext_numperiod)
if len(conditions.pol_ext_spike)>0:
    pol_ext['spike'] = float(conditions.pol_ext_spike)
if len(conditions.pol_ext_rate)>0:
    pol_ext['rate'] = float(conditions.pol_ext_rate)
if len(conditions.pol_ext_start)>0:
    pol_ext['start'] = int(conditions.pol_ext_start)
if len(conditions.pol_ext_species)>0:
    if (conditions.pol_ext_species=='ALL'):
        pol_ext['species'] = ['ALL']
    else:
        pol_ext['species'] = [int(i) for i in conditions.pol_ext_species.split(',')]
        
file_suffix = sgcom.create_file_suffix(algorithm,output_suff + '_' +\
                                       conditions.outsf,year_periods)
fichreport = sgcom.create_fichreport_name(dirsal,filename,file_suffix)

simulation_params = sgcom.SimulationConditions(filename = filename, 
                year_periods = year_periods, 
                hay_foodweb = hay_foodweb, 
                hay_superpredadores = hay_superpredadores,
                data_save = data_save, 
                direntrada = direntrada,
                dirsal = dirsal,
                eliminarenlaces = eliminarenlaces, 
                pl_ext = pl_ext, 
                pol_ext = pol_ext, 
                output_suff = output_suff, 
                fichreport = fichreport, 
                algorithm = algorithm, 
                plants_blossom_prob = plants_blossom_prob, 
                plants_blossom_sd = plants_blossom_sd, 
                plants_blossom_type = plants_blossom_type, 
                blossom_pert_list = blossom_pert_list,
                release = release,
                Bssvar_data = Bssvar_data,
                verbose = verbose)
simulation_params.dirtrabajo = os.getcwd().replace('\\', '/')
simulation_params.exit_on_extinction = conditions.stop_extinction
if verbose:
    print('Report file '+simulation_params.dirtrabajo+'/'+fichreport)
sim_ret_val = b_sim.bino_mutual(sim_cond = simulation_params)
if simulation_params.exit_on_extinction:
    if sim_ret_val.systemextinction:
        print('EXTINCTION')
        sys.exit()
    elif not(verbose):
        print('SURVIVAL')
if conditions.gengraphs:
    sggraph.mutual_render(simulation_params, sim_ret_val, 0,
                                year_periods * sgGL.DAYS_YEAR)
    if hay_foodweb:
                sggraph.food_render(simulation_params, sim_ret_val, 0,
                                year_periods * sgGL.DAYS_YEAR)
    