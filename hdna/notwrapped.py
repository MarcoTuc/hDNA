import os 
import sys 
import gc 
import numpy as np 
import pandas as pd
from hdna import *
from conf import * 
from scipy.stats import gamma

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
        permission = input(f'Folder {RESULTS_DIR} already exists, do you want to overwrite old experiments? [Y,N] ')
        if permission.lower().startswith('y'):
            print('>>>> overwriting old simulations')
            break
        elif permission.lower().startswith('n') or i == 3:
            asknew = input(f'Do you want to make a new folder? [Y,N] ')
            if asknew.lower().startswith('y'):
                control = False
                while True: 
                    NEWEXPNAME = input(f'Write the new experiment name to put inside ./results folder:\n')
                    RESULTS_DIR = f"results/{NEWEXPNAME}"
                    if NEWEXPNAME != EXPNAME: 
                        os.makedirs(RESULTS_DIR)
                        break
                    else: print('Bro, put a new name not the old one... \n')
                break 
            elif permission.lower().startswith('n'): 
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
        sliding_filter=HP['sliding_filter'],
        stacking=OPT['stacking'],
        sliding=HP['sliding_rate'],
        zipping=HP['zipping_rate'],
        celsius=HP['temperature'])

for i, (seq, exp) in enumerate(zip(torun['seq'], torun['expvalue'])):

    print(f'Strand number {i}: {seq}')
    seq = str(seq.strip())      
    exp = float(exp.strip())    

    print(f'Creating network from sequence...')
    A = Strand(MOD, seq)
    B = A.complementary()
    MOD.setgeometry(HP['azimutal_angle'], HP['longitudinal_angle'])

    opts = Options(
        method='direct', 
        runtime=OPT['runtime'], 
        Nsim=OPT['N_simul'], 
        trajstosave=OPT['trajstosave'],
        results_dir=RESULTS_DIR, 
        graphsalone=SOPT['G_saving'],
        stranditer=i)

    print('embedding network into biosimulator network model...')
    simulatore = Simulator(MOD, A, B, options=opts)
    print(simulatore.overview)
    print('start running simulations...')
    results = simulatore.ensemble()
    fpts = simulatore.fpts(results)
    rate = 1/(np.mean(fpts))
    fitgamma = gammafit(fpts)
    fitexp = expfit(fpts)
    histotime(fpts, fitgamma, OPT['runtime'], exp=exp, mod=rate, name=f'{seq}gamma_histo', writepath=simulatore.DIR, theme='dark')
    histotime(fpts, fitexp, OPT['runtime'], exp=exp, mod=rate, name=f'{seq}_exp_histo', writepath=simulatore.DIR, theme='dark')
    percomplot(fpts, writepath=simulatore.DIR, name=f'{seq}_percentplot')
    simulatore.save_graph()
    df = pd.DataFrame.from_dict([simulatore.overview])
    df.drop(['duplex','singlestranded'], axis=1, inplace=True)
    newcols = list(df.columns)
    newvals = list(df.loc[0,df.columns])
    torun.loc[seq, 'computed'] = rate
    torun.loc[seq, newcols] = newvals

    print(f"experimental rate: {'{:e}'.format(exp)}")
    print(f"computed rate:     {'{:e}'.format(rate)}", '\n')


torun.to_csv(f"{RESULTS_DIR}/simulationdata.csv")
valplot(torun, EXPNAME, writepath=RESULTS_DIR, theme='dark')