import networkx as nx
import juliacall 

jl = juliacall.newmodule("hDNA")
jl.seval('using BioSimulator')
jl.seval('using TickTock')

from hdna.kinetwork import Kinetwork
from hdna.model import Model 

class Options(object):
    def __init__(self, method=jl.Direct(), runtime=1, Nmonte=1000):
        self.runtime = runtime
        self.Nmonte = Nmonte
        self.method = method #add some error here if the method is not a julia object or whatever 

class Simulator(object):
    
    def __init__(self, model: Model, kinetwork: Kinetwork, initialamount=2, options=Options):
        
        self.kinet = kinetwork
        self.Graph = kinetwork.Graph
        self.model = model 
        self.options = options 
        self.biosim = jl.Network("biosim")
        self.initialamount = initialamount
        self.add_species()


    def add_species(self):
        """DONE"""
        # self.biosim <= jl.Species(self.tl(str(list(self.Graph.nodes())[0])), self.initialamount)
        ss = self.tl(self.kinet.chamber.singlestranded.structure) 
        self.biosim <= jl.Species(ss, self.initialamount)
        for node in list(self.Graph.nodes())[1:]:
            self.biosim <= jl.Species(self.tl(node))


    def add_reactions(self):
        i = 0  #for now reaction names will just be sequential numbers
        """PSEUDOCODE:
        
        for node in self.Graph:
            get neighborhood of node 
            pairs = node.pairs (how many base pairs the node has)
            state = node.state (what is it)
            for node in neighborhood:
                self.biosim <= jl.Reaction(reaction_name, reaction_rate, reaction_string (use the f"{stuff}" method of Giovanni))
        """
        # start from the neighborhood of the singlestranded state to add nucleations 
        ss = self.kinet.chamber.singlestranded.structure 
        neighbors = nx.neighbors(self.Graph, ss)
        for n in neighbors:
            i += 1                       # TODO RATES
            """ JULIA PATHOLOGY 
                Can't use dots as names for reactions because
                Julia has a pathology with dots and also parenthesis (fuck)
                The most reliable notation is: 

                \   . --> o
                \   ( --> b
                \   ) --> d
                \   + --> A

                """
            ss = self.tl(ss); n  = self.tl(n)
            name = f"f{i}"; rule = f"{ss} + {ss} --> {n}"
            self.biosim <= jl.Reaction(name, 0, rule)
            name = f"b{i}"; rule = f"{n} --> {ss} + {ss}"
            self.biosim <= jl.Reaction(name, 0, rule)

        # now make a graph without the singlestranded state for all subsequent transitions 
        subgraph = self.Graph.copy().remove_node(ss)
        for n1, n2, data in list(subgraph.edges.data()):
            i += 1                       # TODO RATES (data will be used to get rates etc)
            n1 = self.tl(n1)
            n2 = self.tl(n2)
            name = f'f{i}'; rule = f'{n1} --> {n2}'
            self.biosim <= jl.Reaction(name, 0, rule)
            name = f'b{i}'; rule = f'{n2} --> {n1}'
            self.biosim <= jl.Reaction(name, 0, rule)




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
    
    def trajectory(self):
        pass
        """ Create a method for appending a simulation trajectory to each simulation result """

    ######################
    ### Helper Methods ###
    ######################

    def tl(self, s):
        """ Used to fix julia pathology """
        table = str.maketrans({'.':'o', '(':'b', ')':'d', '+':'A'})
        return s.translate(table)

    # DEPRECATED
    # def nucleation_filter(self):
    #     return nx.edge_subgraph(self.Graph, edges=
    #     [e for e, attrdict in self.Graph.edges.items() if attrdict['kind' == 'on_nucleation'] or attrdict['kind' == 'off_nucleation']])