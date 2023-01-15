import os
os.environ['PYTHON_JULIACALL_THREADS'] = '6'
os.environ['PYTHON_JULIACALL_PROCS'] = '6'

import networkx as nx

import juliapkg
juliapkg.resolve()

import juliacall 
jl = juliacall.Main

jl.seval('using BioSimulator')
jl.seval('using TickTock')

print(jl.CONFIG)
print(jl.Threads.nthreads())

from hdna.kinetwork import Kinetwork, Kinetics
from hdna.model import Model, Geometry

class Options(object):
    def __init__(self, method="direct", runtime=1e-5, Nmonte=500):
        self.runtime = runtime
        self.Nmonte = Nmonte

        methods = {"direct": jl.Direct()}

        self.method = methods[method] #add some error here if the method is not a julia object or whatever 

class Simulator(object):
    
    def __init__(self, model: Model, kinetwork: Kinetwork, kinetics: Kinetics, initialamount=2, options=Options()):
        
        self.model = model 
        self.kinet = kinetwork
        self.Graph = kinetwork.Graph

        self.kinetics = kinetics 

        self.options = options 
        self.biosim = jl.Network("biosim")

        self.initialamount = initialamount
        self.add_species()
        self.add_reactions()


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
        SS = self.kinet.chamber.singlestranded.structure 
        ss = self.tl(SS)
        neighbors = list(nx.neighbors(self.Graph, SS))
        sliding_factor = 1/len(neighbors)
        for n in neighbors:
            i += 1
            """ JULIA PATHOLOGY 
                Can't use dots as names for reactions because
                Julia has a pathology with dots and also 
                parenthesis (fuck) the most reliable notation is: 
                |   . --> o     |
                |   ( --> b     |
                |   ) --> d     |
                |   + --> A     |                   """
            neigh  = self.tl(n)

            # FORWARD BIMOLECULAR NUCLEATION
            name = f"f_nucleation{i}"; rule = f"{ss} + {ss} --> {neigh}"
            kf = self.kinetics.georate * sliding_factor
            self.biosim <= jl.Reaction(name, kf, rule)
            self.Graph.edges[SS,n]['kf'] = kf   # Also add the rate as a new property to the graph, it will be useful for visualization

            # BACKWARD UNIMOLECULAR DISSOCIATION            
            name = f"b_nucleation{i}"; rule = f"{neigh} --> {ss} + {ss}"
            DG = self.Graph.nodes[n]['object'].G
            kb = self.kinetics.k_back(kf, DG)
            self.biosim <= jl.Reaction(name, kb, rule)
            self.Graph.edges[SS,n]['kb'] = kb   # Also add the rate as a new property to the graph, it will be useful for visualization

        # now make a graph without the singlestranded state for all subsequent transitions 
        subgraph = self.Graph.copy()
        subgraph.remove_node(SS)
        
        for n1, n2, data in list(subgraph.edges.data()):
            i += 1                       # TODO RATES (data will be used to get rates etc)
            
            if data['kind'] in ['sliding', 'sliding-end']:
                l1 = self.Graph.nodes[n1]['object'].total_nucleations
                l2 = self.Graph.nodes[n2]['object'].total_nucleations
                if l1 < l2:             # Left slidings
                    n1o = self.tl(n1)
                    n2o = self.tl(n2)
                else:                   # Right slidings
                    n1o = self.tl(n2)
                    n2o = self.tl(n1)
            else:
                n1o = self.tl(n1)
                n2o = self.tl(n2)

            # FORWARD SLIDINGS AND ZIPPINGS
            name = f"f_{data['kind']}_{i}"; rule = f"{n1o} --> {n2o}"
            kf = self.kinetics.generalforward(data['kind'])
            self.biosim <= jl.Reaction(name, kf, rule)
            self.Graph.edges[n1,n2][f'kf'] = kf   # Also add the rate as a new property to the graph, it will be useful for visualization

            # BACKWARD SLIDINGS AND ZIPPINGS 
            name = f"b_{data['kind']}_{i}"; rule = f'{n2o} --> {n1o}'
            DG = self.Graph.nodes[n2]['object'].G - self.Graph.nodes[n1]['object'].G
            if self.Graph.nodes[n1]['object'].total_nucleations > self.Graph.nodes[n2]['object'].total_nucleations: DG = - DG
            kb = self.kinetics.k_back(kf, DG)
            self.biosim <= jl.Reaction(name, kb, rule)
            self.Graph.edges[n1,n2][f'kb'] = kb   # Also add the rate as a new property to the graph, it will be useful for visualization

            """ Note that for slidings coming from the right like: .......((((+.......)))) <--> ........(((+........))) 
                the forward rate will always denote the passage from the strand with less basepairs with the one with more.
                By now the network visualization is kinda broken but who cares, rates are where they should be. 
                This happens because of how nodes are ordered inside kinetwork's Graph """

    def simulation(self):
        """ run the simulation and return the results """
        sim = jl.simulate(self.biosim, jl.NextReaction(), tfinal = self.options.runtime)
        return {'c':jl.hcat(*sim.u).__array__(), 't':sim.t}
    
    def ensemble(self):
        state, model = jl.parse_model(self.biosim)
        return [jl.simulate(state, model, self.options.method, tfinal = self.options.runtime) for _ in range(self.options.Nmonte)]

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
        """ Create a method for appending a simulation trajectory to each simulation result """
        pass

    ######################
    ### Helper Methods ###
    ######################

    def tl(self, s):
        """ Used to fix julia pathology """
        table = str.maketrans({'.':'o', '(':'b', ')':'d', '+':'A'})
        return s.translate(table)
    
    def lt(self,s):
        table = str.maketrans({'o':'.', 'b':'(', 'd':')', 'A':'+'})
        return s.translate(table)



    # DEPRECATED
    # def nucleation_filter(self):
    #     return nx.edge_subgraph(self.Graph, edges=
    #     [e for e, attrdict in self.Graph.edges.items() if attrdict['kind' == 'on_nucleation'] or attrdict['kind' == 'off_nucleation']])