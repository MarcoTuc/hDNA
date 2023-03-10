import os 
import sys 
import numpy as np 
import pandas as pd

EXPNAME = 'goodmorning-5'
RESULTS_DIR = f"results/{EXPNAME}"

HP = {
    #model free parameters  
    'minimum_nucleation': 1,
    'sliding_cutoff':     200,
    'sliding_filter':     4,
    'zipping_rate':       8e7,
    'sliding_rate':       2e7,
    #temperature
    'temperature':        25,       #### HERTEL EXPERIMENTAL TEMPERATURE 
    #angles
    'azimutal_angle':     120,
    'longitudinal_angle': 270,
}

OPT = {
    #simulation options
    'runtime': 4e-6,
    'N_simul': 1000,
    'trajstosave': 35,
    #nupack options
    'stacking': 'nostacking'
}

SOPT = {    
    #datasaving options 
    'G_saving': 'strand_folder'
}
