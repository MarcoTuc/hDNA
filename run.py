import numpy as np 
import pandas as pd
from hdna import *


# Import experimental data from Hertel 
expdata = pd.read_csv('./data/herteldata.csv', names=['seq', 'expvalue'])
# Clean the dataframe 
expdata = expdata.drop(0)
expdata['expvalue'] = ['{:e}'.format(float(e)) for e in expdata['expvalue']]

# Actual computation 

rates = []

model = Model('dna', '3D')

for i, (seq, exp) in enumerate(zip(expdata['seq'], expdata['expvalue'])):
    seq = str(seq.strip())      # make sure they are alright
    exp = float(exp.strip())    #
    print(f'Creating network from string {seq}...')
    A = Strand(model, seq)
    B = A.complementary()
    kinet = Kinetwork(model, A, B, 1)
    geo = Geometry(120, 360)
    K = Kinetics(model, kinet, geo)
    opts = Options(method='direct', runtime=3e-6, Nmonte=20000)
    print('embedding strand into biosimulator network model...')
    simulatore = Simulator(model, kinet, K, options=opts)
    print('start running simulations...')
    results = simulatore.ensemble()
    mfpt = simulatore.mfpts(results)
    rates.append(1/mfpt)
    print(f"experimental rate: {'{:e}'.format(exp)}")
    print(f"computed rate:     {'{:e}'.format(1/mfpt)}")
    del results 