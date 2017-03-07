''' This script executes a simulation, takes the final population values
and creates new input matrixes with those values as initial ones

Example: sigmund_createinput_atmax.py -fichout output/output_data_M1_A1_all_Verhulst_a_populations__100.txt '''

import argparse
import os
import sigmund_GLOBALS as sgGL
import sigmund_common as sgcom
import sgtools_common as sgtcom
import math

def rewrite_input_file(listmax, inputfile, listaent):
    listaent[-5] = listmax[:]
    cuts = inputfile.split('_')
    filemax = '_'.join(cuts[0:-1]) + '_MAX_' + cuts[-1]
    print('Creating '+sgGL.INPUTFILES_PATH+str(filemax))
    salida = open(sgGL.INPUTFILES_PATH + filemax, 'w', encoding='utf-8')
    for linea in listaent:
        textosal = '\t'.join(linea)
        textosal = textosal + '\n'
        salida.write(textosal)
    salida.close()

def read_results_file(filename):
    with open(filename) as myfile:
        ultlin = list(myfile)[-1]
    ultlin = ultlin.split(',')
    listmax = [(i.split('\n')[0]).split('.')[0] for i in ultlin]
    return listmax


parser = argparse.ArgumentParser()
parser.add_argument('-fichout', action='store', dest='fichout_name',default = '',
                    help='Populations output file with _a suffix')
                        
conditions = parser.parse_args()

filename = conditions.fichout_name
filename_a = filename
print('filename '+str(filename_a))
listmax_a = read_results_file(filename_a)
print('Plants maxima '+str(listmax_a))

filename_b = filename.replace('_a_populations', '_b_populations')
listmax_b = read_results_file(filename_b)
print('Pollinators maxima '+str(listmax_b))

inputname = filename.replace('\\', '/ ').split('/')[-1]
inputname = inputname.split('Verhulst')[0]
inputname = inputname.split('output_data_')[-1]
inputfile_a = inputname+'a.txt'
inputfile_b = inputname+'b.txt'
dirtra = os.getcwd().replace('\\', '/')

''' Read input matrixes '''
listaent_a = sgcom.dlmreadlike(inputfile_a,sgGL.INPUTFILES_PATH)
pop_ini_a = listaent_a[-5][:]
listaent_b = sgcom.dlmreadlike(inputfile_b,sgGL.INPUTFILES_PATH)
pop_ini_b = listaent_b[-5][:]

rewrite_input_file(listmax_a, inputfile_a, listaent_a)
rewrite_input_file(listmax_b, inputfile_b, listaent_b)