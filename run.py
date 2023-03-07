
import pandas as pd 
import numpy as np 

from hdna import *
from tqdm import tqdm

from scipy.optimize import dual_annealing

EXPNAME = 'sumofinchandpk17'

notes = """

manually lowered fwd for slidings with only two base pairs by a factor of 100

"""

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
MOD.setgeometry(theta=120, phi = 270)
 
MOD.alpha = 1
MOD.gamma = 0
MOD.kappa = 1

OPT = Options(Nsim=5000)

H = HDNA(data, EXPNAME, model=MOD, options=OPT)

with open(f'results/{EXPNAME}/notes.txt', 'w') as savenote:
    savenote.write(notes)
    savenote.close()


zipping = 10e7
sliding = 2e5 #zipping*np.exp(-(2)/(CONST.R*MOD.kelvin))

H.run([zipping, sliding])