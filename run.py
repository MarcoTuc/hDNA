import numpy as np 
import pandas as pd
from hdna import *

RESULTS_DIR = "results/min4long"

# Import experimental data from Hertel 
expdata = pd.read_csv('./data/herteldata.csv', names=['seq', 'expvalue'])
# Clean the dataframe 
expdata = expdata.drop(0)
expdata['expvalue'] = ['{:e}'.format(float(e)) for e in expdata['expvalue']]

# Actual computation 

rates = []

model = Model('dna', '3D')

for i, (seq, exp) in enumerate(zip(expdata['seq'], expdata['expvalue'])):
    print(f'Strand number {i+1}: {seq}')
    seq = str(seq.strip())      # make sure they are alright
    exp = float(exp.strip())    #
    print(f'Creating network from sequence...')
    A = Strand(model, seq)
    B = A.complementary()
    kinet = Kinetwork(model, A, B, 4)
    geo = Geometry(120, 360)
    K = Kinetics(model, kinet, geo)
    opts = Options(method='direct', runtime=3e-6, Nmonte=2000)
    print('embedding network into biosimulator network model...')
    simulatore = Simulator(model, kinet, K, options=opts)
    print('start running simulations...')
    results = simulatore.ensemble(write_csv=True, folder=RESULTS_DIR, stranditer = i+1)
    mfpt = simulatore.mfpts(results)
    rates.append(1/mfpt)
    print(f"experimental rate: {'{:e}'.format(exp)}")
    print(f"computed rate:     {'{:e}'.format(1/mfpt)}")
    del results 

expdata['computedrates'] = rates
expdata.to_csv(RESULTS_DIR)

#TODO
""" Create a routine for a more easy identification of strands
    (maybe add a number to strand folders so to find them easy)
    Also put in each strand folder an html file containing its 
    kinetwork for visualization and error checking. 
    Also for error checking put some csv with nodes and edges of the 
    network and with the same information but from the biosim model. 
"""