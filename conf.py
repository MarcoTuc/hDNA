import os 
import sys 
import numpy as np 
import pandas as pd
from hdna.logger import Tee 

EXPNAME = 'fullgraph_0'
RESULTS_DIR = f"results/{EXPNAME}"

HP = {
    #model free parameters  
    'minimum_nucleation': 1,
    'sliding_cutoff':     2,
    'zipping_rate':       2e8,
    'sliding_rate':       5e5,
    #temperature
    'temperature':        25,       #### HERTEL EXPERIMENTAL TEMPERATURE 
    #angles
    'azimutal_angle':     120,
    'longitudinal_angle': 270,
}

OPT = {

    #simulation options
    'runtime': 5e-6,
    'N_simul': 35,
    'trajstosave': 35,

    #nupack options
    'stacking': 'nostacking'
}

SOPT = {    
    
    #datasaving options 
    'G_saving': 'strand_folder'
}
