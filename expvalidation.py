import os 
os.environ['JULIA_NUM_THREADS'] = '1'
os.environ['PYTHON_JULIACALL_THREADS'] = '1'
os.environ['PYTHON_JULIACALL_PROCS'] = '1'

import juliacall 

# import pandas as pd
# from hdna import *

# expdata = pd.read_csv('./data/herteldata.csv')

# model = Model('dna', '3D')

# data = expdata.iloc[0]

# print(f'Simulating string {(data[0]).strip()}')
# A = Strand(model, str(data[0]).strip())
# B = A.complementary()
# kinet = Kinetwork(model, A, B, 3)
# geo = Geometry(180, 360)
# K = Kinetics(model, kinet, geo)
# opts = Options(method='direct', runtime=2e-6, Nmonte=10000)
# simulatore = Simulator(model, kinet, K, options=opts)

# sim = simulatore.ensemble()
