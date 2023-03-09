
import pandas as pd 
import numpy as np 

from hdna import *
from tqdm import tqdm

from scipy.optimize import dual_annealing

EXPNAME = 'HATA1'

notes = """

HATA DATA 

"""

# Import experimental data from Hertel 
expdata = pd.read_csv('./data/hatadata.csv')
# Clean the dataframe 
expdata = expdata[['rate 10^5 /M /s', 'sequence1']]
expdata = expdata.rename(columns={'rate 10^5 /M /s':'experimental','sequence1': 'sequences'})
expdata['experimental'] = ['{:3e}'.format(float(e*1e5)) for e in expdata['experimental']]

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

OPT = Options(Nsim=5000)

H = HDNA(data, EXPNAME, model=MOD, options=OPT)

with open(f'results/{EXPNAME}/notes.txt', 'w') as savenote:
    savenote.write(notes)
    savenote.close()


zipping = 10e7
sliding = 1.5e5 #zipping*np.exp(-(2)/(CONST.R*MOD.kelvin))

H.run([zipping, sliding])