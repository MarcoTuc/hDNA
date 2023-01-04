import networkx as nx

from classes.chamber import Chamber
from classes.model import Model
 
class Kinnet(object):
    """Given a chamber it will generate the corresponding
    Kinetic network to be later passed to the simulator"""

    def __init__(self, model: Model, chamber: Chamber, minimum_nucleation):

        self.model = model 
        self.chamber = chamber

        """ offnodes are chamber's complexes """
        self.offnodes = chamber.offcores

        """ onnodes are chamber's complexes and each of these 
            complexes has in itself a zipping trajectory which is 
            a property of the complex and not of the chamber """
        self.onnodes = chamber.oncores

        
