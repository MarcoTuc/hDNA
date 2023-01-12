import juliacall 

jl = juliacall.newmodule("hDNA")
jl.seval('using BioSimulator')
jl.seval('using TickTock')

from hdna.kinetwork import Kinetwork
from hdna.model import Model 

class Options(object):
    def __init__(self, method, runtime, Nmonte):
        self.runtime = runtime
        self.Nmonte = Nmonte
        self.method = method #add some error here if the method is not a julia object or whatever 

class Simulator(object):
    
    def __init__(self, model: Model, kinetwork: Kinetwork, initialamount, options):
        
        self.Graph = kinetwork.Graph
        self.model = model 
        self.options = options 
        self.biosim = jl.Network("biosim")
        self.initialamount = initialamount


    def add_species(self):
        """PSEUDOCODE
        for node in self.Graph:
            self.biosim <= jl.Species(node.structure, if it is singlestranded add the self.initial)    
        """
        


    def add_reactions(self):
        """PSEUDOCODE:
        
        for node in self.Graph:
            get neighborhood of node 
            pairs = node.pairs (how many base pairs the node has)
            state = node.state (what is it)
            for node in neighborhood:
                self.biosim <= jl.Reaction(reaction_name, reaction_rate, reaction_string (use the f"{stuff}" method of Giovanni))
        """
        pass
    
    def run_simulation(self):
        """
        return [jl.simulate(self.biosim, self.options.method, tfinal=self.options.runtime) for i in range(self.options.Nmonte)]
        """
        pass
    
    def mfpts(self, simresults):
        """
        node referring to Duplex should be the last one so the result spot for duplex should just be something like simresults[-1] i hope 
        first i need to understand how the simresult looks like out of the juliacalled biosimlator 

        PSEUDOCODE:

        def condition(x):
            return 1 == y (when the value of duplex gets to one for the firs time we get the mfpt)
        def findmfpt(sim):
            index_r = findfirst(condition, result[:,duplexindex])
            firstpassage = etc etc 

        mfpts = []
        time = simresults.time
        for sim in simresults:
            mfpts.append(self.findmfpt(sim)) make another method that just gets the mfpt for a single simulation result 
        """
    
