import networkx as nx 
from .chamber import Chamber
from .model import Model 
 
class Kinnet(object):
    """Given a chamber it will generate the corresponding
    Kinetic network to be later passed to the simulator"""

    def __init__(self, model: Model, chamber: Chamber, minimum_nucleation):

        self.model = model 
        self.chamber = chamber

        self.nodes = self.chamber.compute_cores(minimum_nucleation)