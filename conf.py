import os 
import sys 
import numpy as np 
import pandas as pd
from hdna.logger import Tee 

EXPNAME = 'nosliding_1'
RESULTS_DIR = f"results/{EXPNAME}"

HP = {
    
    #model free parameters  
    'minimum_nucleation': 4,
    'sliding_cutoff':     0,
    'zipping_rate':       2e9,
    'sliding_rate':       2e7,
    
    #temperature
    'temperature':        25,       #### HERTEL EXPERIMENTAL TEMPERATURE 
    
    #angles
    'azimutal_angle':     120,
    'longitudinal_angle': 270,
}

OPT = {

    #simulation options
    'runtime': 5e-6,
    'N_simul': 1500,
    'trajstosave': 35,
}

SOPT = {    
    
    #datasaving options 
    'G_saving': 'strand_folder'
}
