import networkx as nx

import os
os.environ['PYTHON_JULIACALL_THREADS'] = '4'
os.environ['PYTHON_JULIACALL_PROCS'] = '4'

# import juliapkg
# juliapkg.resolve()

import juliacall 


# jl.seval('@everywhere using BioSimulator')
# jl.seval('@everywhere using TickTock')

# print(juliacall.CONFIG)
# print(jl.Threads.nthreads())