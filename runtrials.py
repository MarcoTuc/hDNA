from hdna.strand import Strand
from hdna.model import Model, Geometry
from hdna.kinetwork import Kinetwork, Kinetics
from hdna.simulator import Simulator

print('imports done')

model = Model('dna', '3D')
mincore = 3

a = Strand(model, 'AAAAAAAAAAA')
b = Strand(model, 'TTTTTTTTTTT')

kgraph = Kinetwork(model, a, b, mincore)
geometry = Geometry(180, 180)
kinetics = Kinetics(model, kgraph, geometry)
simulator = Simulator(model, kgraph, kinetics)

print('simulator made')

simulator.simulation()