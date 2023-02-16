
import pandas as pd 

from hdnawrapper import HDNA
from tqdm import tqdm

from scipy.optimize import dual_annealing

EXPNAME = 'firstrial'

# Import experimental data from Hertel 
expdata = pd.read_csv('./data/herteldata.csv', names=['sequences', 'experimental'])
# Clean the dataframe 
expdata = expdata.drop(0)
expdata['experimental'] = ['{:e}'.format(float(e)) for e in expdata['experimental']]

limit = 2
data = expdata.copy().iloc[:limit]
data['index'] = data.index 
data.set_index(data['sequences'], inplace=True)

H = HDNA(data, EXPNAME)
results = dual_annealing(H.run, [2e7, 2e7], niter=5, stepsize=2e7, T=3e15)