import os 
import sys 
import numpy as np 
import pandas as pd
from hdna import *

# Import experimental data from Hertel 
expdata = pd.read_csv('./data/herteldata.csv', names=['seq', 'expvalue'])
# Clean the dataframe 
expdata = expdata.drop(0)
expdata['expvalue'] = ['{:e}'.format(float(e)) for e in expdata['expvalue']]

limit = len(expdata)
torun = expdata.copy().iloc[:limit]

EXPNAME = 'MN5Z9S5_(180,360)'
RESULTS_DIR = f"results/{EXPNAME}"

if os.path.isdir(RESULTS_DIR): 
    i = 0
    while True: 
        i += 1
        permission = input('Folder already exists, do you want to overwrite old experiments? [Y,N]')
        if permission.lower().startswith('y'):
            print('>>>> overwriting old simulations')
            break
        elif permission.lower().startswith('n') or i == 3:
            print(">>>> stopping the program")
            sys.exit()
        print("yes or not?") 
else:
    os.makedirs(RESULTS_DIR)

HP = {
    
    #model free parameters  
    'minimum_nucleation': 5,
    'zipping_rate':       2e9,
    'sliding_rate':       5e5,
    
    #temperature
    'temperature':        37,
    
    #angles
    'azimutal_angle':     180,
    'longitudinal_angle': 360,
}

OPT = {

    #simulation options
    'runtime': 5e-6,
    'N_simul': 3000
}

hyperparams = pd.DataFrame.from_dict([dict(**HP,**OPT)]).T
hyperparams.rename(columns={np.int64(0):'values'}, inplace=True)
hyperparams.index.rename('hyperparameters', inplace=True)
hyperparams.to_csv(f'{RESULTS_DIR}/hyperparameters.csv')

# Actual computation 
rates = []
model = Model('dna', '3D', celsius=HP['temperature'])

for i, (seq, exp) in enumerate(zip(torun['seq'], torun['expvalue'])):
    print(f'Strand number {i+1}: {seq}')
    seq = str(seq.strip())      # make sure they are alright
    exp = float(exp.strip())    #
    print(f'Creating network from sequence...')
    A = Strand(model, seq)
    B = A.complementary()
    kinet = Kinetwork(model, A, B, HP['minimum_nucleation'])
    geo = Geometry(HP['azimutal_angle'], HP['longitudinal_angle'])
    K = Kinetics(model, kinet, geo)
    K.set_slidingrate(HP['sliding_rate'])
    K.set_zippingrate(HP['zipping_rate'])
    opts = Options(method='direct', runtime=OPT['runtime'], Nsim=OPT['N_simul'], results_dir=RESULTS_DIR, stranditer=i)
    print('embedding network into biosimulator network model...')
    simulatore = Simulator(model, kinet, K, options=opts)
    print('start running simulations...')
    results = simulatore.ensemble()
    mfpt = simulatore.mfpts(results)
    rates.append(1/mfpt)
    print(f"experimental rate: {'{:e}'.format(exp)}")
    print(f"computed rate:     {'{:e}'.format(1/mfpt)}", '\n')
    del results

torun['computed'] = rates
torun.to_csv(f"{RESULTS_DIR}/simulationdata.csv")
valplot(torun, EXPNAME, writepath=RESULTS_DIR, theme='dark')


#TODO
""" Also for error checking put some csv with nodes and edges of the 
    network and with the same information but from the biosim model. 
"""