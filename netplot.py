import matplotlib.pyplot as plt
import hvplot.networkx as hvnx
from pyvis.network import Network
import nupack as nu
import networkx as nx 
import pandas as pd 
from classes.complex import Complex
from classes.strand import Strand
from classes.model import Model
from classes.chamber import Chamber
from classes.kinetwork import Kinetwork

model = Model('dna', '3D')
mincore = 4

a = Strand(model, 'AAAAAAAAAAAAAAAAAA')
b = Strand(model, 'TTTTTTTTTTTTTTTTTT')

kgraph = Kinetwork(model, a, b, mincore)

G = kgraph.Graph

net = Network(notebook=True,cdn_resources='remote')
net.from_nx(G)
net.show('G.html')