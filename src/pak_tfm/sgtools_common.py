#!/usr/bin/env python3
import pickle
import sigmund_GLOBALS as sgGL
import sigmund_common as sgcom

def read_file_input(simulation_file):
    path_sims = sgGL.INPUTFILES_PATH+sgGL.SIMFILES_PATH
    try:
        fh = open(path_sims+simulation_file, "rb")
    except IOError:
        print(simulation_file+' not found')
        return
    else:
        fh.close()
        with open(path_sims+simulation_file, 'rb') as f:
            storedsim = pickle.load(f)
            return(storedsim)