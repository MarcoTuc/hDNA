import os 
import sys 
import gc 
import numpy as np 
import pandas as pd
from hdna import *
from conf import * 

# Import experimental data from Hertel 
expdata = pd.read_csv('./data/herteldata.csv', names=['seq', 'expvalue'])
# Clean the dataframe 
expdata = expdata.drop(0)
expdata['expvalue'] = ['{:e}'.format(float(e)) for e in expdata['expvalue']]

limit = len(expdata)
torun = expdata.copy().iloc[:limit]
torun['index'] = torun.index 
torun.set_index(torun['seq'], inplace=True)

# Directory Check 
if os.path.isdir(RESULTS_DIR): 
    i = 0
    while True: 
        i += 1
        permission = input(f'Folder {RESULTS_DIR} already exists, do you want to overwrite old experiments? [Y,N]')
        if permission.lower().startswith('y'):
            print('>>>> overwriting old simulations')
            break
        elif permission.lower().startswith('n') or i == 3:
            print(">>>> stopping the program")
            sys.exit()
        print("yes or not?") 
else:
    os.makedirs(RESULTS_DIR)

# Leave a csv with hyperparameters 
hyperparams = pd.DataFrame.from_dict([dict(**HP,**OPT)]).T
hyperparams.rename(columns={np.int64(0):'values'}, inplace=True)
hyperparams.index.rename('hyperparameters', inplace=True)
hyperparams.to_csv(f'{RESULTS_DIR}/hyperparameters.csv')

# export console output to a txt 
f = open(f'{RESULTS_DIR}/console.txt', 'w')
sys.stdout = Tee(sys.stdout, f)

# Actual computation 
MOD = Model('dna', '3D', 
        min_nucleation=HP['minimum_nucleation'], 
        sliding_cutoff=HP['sliding_cutoff'],
        celsius=HP['temperature'])

for i, (seq, exp) in enumerate(zip(torun['seq'], torun['expvalue'])):

    print(f'Strand number {i}: {seq}')
    seq = str(seq.strip())      
    exp = float(exp.strip())    

    print(f'Creating network from sequence...')
    A = Strand(MOD, seq)
    B = A.complementary()
    kinet = Kinetwork(MOD, A, B)
    geo = Geometry(HP['azimutal_angle'], HP['longitudinal_angle'])
    K = Kinetics(MOD, kinet, geo)
    K.set_slidingrate(HP['sliding_rate'])
    K.set_zippingrate(HP['zipping_rate'])
    print(kinet.overview)

    opts = Options(
        method='direct', 
        runtime=OPT['runtime'], 
        Nsim=OPT['N_simul'], 
        trajstosave=OPT['trajstosave'],
        results_dir=RESULTS_DIR, 
        graphsalone=SOPT['G_saving'],
        stranditer=i)

    print('embedding network into biosimulator network model...')
    simulatore = Simulator(MOD, kinet, K, options=opts)
    print('start running simulations...')
    results = simulatore.ensemble()
    mfpt = simulatore.mfpts(results)

    df = pd.DataFrame.from_dict([simulatore.overview])
    df.drop(['duplex','singlestranded'], axis=1, inplace=True)
    newcols = list(df.columns)
    newvals = list(df.loc[0,df.columns])
    torun.loc[seq, 'computed'] = 1/mfpt
    torun.loc[seq, newcols] = newvals

    print(f"experimental rate: {'{:e}'.format(exp)}")
    print(f"computed rate:     {'{:e}'.format(1/mfpt)}", '\n')


torun.to_csv(f"{RESULTS_DIR}/simulationdata.csv")
valplot(torun, EXPNAME, writepath=RESULTS_DIR, theme='dark')


#TODO
""" Also for error checking put some csv with nodes and edges of the 
    network and with the same information but from the biosim model. 
"""