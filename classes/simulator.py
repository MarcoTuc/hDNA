import juliacall 

jl = juliacall.newmodule("hDNA")
jl.seval('using BioSimulator')

from classes.kinetwork import Kinetwork

class Simulator(object):
    
    def __init__(self, kinetwork: Kinetwork):
        
        self.kinetwork = kinetwork
