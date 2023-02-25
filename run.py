
import pandas as pd 

from hdna import *
from tqdm import tqdm

from scipy.optimize import dual_annealing

EXPNAME = 'lengthdependant'

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
    min_nucleation=1,
    normalizeback=False)
MOD.setgeometry(theta=120, phi = 270)
 
MOD.alpha = 0
MOD.gamma = 0
MOD.kappa = 1/6

OPT = Options(Nsim=2500)

H = HDNA(data, EXPNAME, model=MOD, options=OPT)
# bounds = [(2e7, 2e8), (2e6, 2e8)]
# results = dual_annealing(H.run, bounds, maxiter=5, initial_temp=500)

zipping = 8.5e7
sliding = 5e6

H.run([zipping, sliding])