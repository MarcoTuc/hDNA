
import pandas as pd 

from hdna import *
from tqdm import tqdm

from scipy.optimize import dual_annealing

EXPNAME = '1e2'

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
    space_dimensionality='2D',
    stacking='nostacking',
    min_nucleation=1)
MOD.setgeometry(theta=180, phi = 270)
 
MOD.alpha = 1
MOD.gamma = 0
MOD.kappa = 1

OPT = Options(Nsim=200,
              runtime=1e-5)

H = HDNA(data, EXPNAME, model=MOD, options=OPT)
# bounds = [(2e7, 2e8), (2e6, 2e8)]
# results = dual_annealing(H.run, bounds, maxiter=5, initial_temp=500)

zipping = 7.5e7
sliding = 1.5e5

H.run([zipping, sliding])