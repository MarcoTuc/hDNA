
import pandas as pd 
import numpy as np 

from hdna import *
from tqdm import tqdm


from scipy.optimize import dual_annealing

EXPNAME = 'inchpseudo9'

notes = """
let's see"""

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
MOD.setgeometry(theta=90, phi = 120)
 
MOD.alpha = 1
MOD.gamma = 0
MOD.kappa = 1

OPT = Options(Nsim=2500)

H = HDNA(data, EXPNAME, model=MOD, options=OPT)
# bounds = [(2e7, 2e8), (2e6, 2e8)]
# results = dual_annealing(H.run, bounds, maxiter=5, initial_temp=500)
with open(f'results/{EXPNAME}/notes.txt', 'w') as savenote:
    savenote.write(notes)
    savenote.close()

zipping = 7.5e7
sliding = zipping*np.exp(-(1.5)/(CONST.R*MOD.kelvin))

H.run([zipping, sliding])

