
import pandas as pd 
import os 
from hdna import *
from scipy.optimize import dual_annealing

from multiprocessing import Process, Pool

NUM_PROCESSES = os.cpu_count()

def hdnarun(expname, zipping, sliding):

    # Import experimental data from Hertel 
    expdata = pd.read_csv('./data/herteldata.csv', names=['sequences', 'experimental'])
    # Clean the dataframe 
    expdata = expdata.drop(0)
    expdata['experimental'] = ['{:e}'.format(float(e)) for e in expdata['experimental']]

    limit = len(expdata)
    data = expdata.copy().iloc[:limit]
    data['index'] = data.index 
    data.set_index(data['sequences'], inplace=True)

    MOD = Model(
        stacking='nostacking',
        min_nucleation=1)
    MOD.setgeometry(theta=180, phi = 270)
    
    MOD.alpha = 1
    MOD.gamma = 0
    MOD.kappa = 1

    OPT = Options(Nsim=5000)

    H = HDNA(data, expname=expname, model=MOD, options=OPT)

    H.run([zipping, sliding])



argz = (['exp1', 7.5e7, 1.5e7], ['exp2', 7e7, 1.5e7], ['exp3', 6e7, 2e7])

if __name__== "__main__":
    with Pool(processes=NUM_PROCESSES) as pool:
        pool.starmap_async(hdnarun, argz)
