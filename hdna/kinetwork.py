import networkx as nx
import numpy as np 
import pandas as pd 

from itertools import pairwise, combinations, tee, permutations

from .complex import Complex
from .model import Model
from .strand import Strand, Structure
from .kinetics import Kinetics
 
class Kinetwork(object):

    """
    Given a chamber it will generate the corresponding
    Kinetic network to be later passed to the simulator.

    The Zipping Graph is the composition of single zipping graphs ---> See zippingraph() method
    corresponding to each oncore taken from the chamber.
    The Sliding Graph is the composition of single sliding graphs ---> See slidingraph() method
    """


##########################################################################

    def __init__(self, model: Model, s1: Strand, s2: Strand, clean=False):

        self.model = model 
        self.s1 = s1        #53
        self.s2 = s2.invert #35
        self.min_nucleation = self.model.min_nucleation
        self.sliding_cutoff = self.model.sliding_cutoff
        
        self.kinetics = Kinetics(self.model, s1, s2)
        self.kinetics.set_slidingrate(self.model.sliding)
        self.kinetics.set_zippingrate(self.model.zipping)
        
        self.ssobj = Complex(self.model, self.s1, self.s2, state='singlestranded', structure=self.simplex, dpxdist=max([self.s1.length, self.s2.length]))
        self.dxobj = Complex(self.model, self.s1, self.s2, state='duplex', structure=self.duplex, dpxdist=0)

        if self.model.space_dimensionality == '3D':
            self.nmethod = self.kinetics.nucleation
            self.zmethod = self.kinetics.metropolis
            self.smethod = self.kinetics.kawasaki
        else: 
            self.nmethod = self.kinetics.topbotnucleation
            self.zmethod = self.kinetics.metropolis
            self.smethod = self.kinetics.kawasaki
            
        if not clean: 

            self.completegraph()

            #get an overview of the network 
            self.possible_states = [ None,
                                    'singlestranded',
                                    'duplex',
                                    'zipping',
                                    'on_nucleation',
                                    'off_nucleation', 
                                    'backfray',
                                    'sliding' ]   

            countslist = [list(dict(self.DG.nodes.data('state')).values()).count(state) for state in self.possible_states]
            self.overview = dict(zip(self.possible_states,countslist))


##########################################################################
##########################################################################

    def completegraph(self):

        self.DG = nx.DiGraph()

        self.DG.add_node(self.simplex, 
                         obj = self.ssobj,
                         pairs = 0, 
                         state = 'singlestranded',
                         dpxdist = self.ssobj.dpxdist,
                         fre = 0)
        self.DG.add_node(self.duplex, 
                         obj = self.dxobj, 
                         pairs = self.dxobj.total_nucleations, 
                         state = 'duplex',
                         dpxdist = 0,
                         fre = self.dxobj.G)
        
        self.get_graph()
        self.connect_slidings()


##########################################################################
##########################################################################

    def get_graph(self, verbose=False):
        self.sldbranches = []
        ss = self.simplex.split('+')[0] 
        num = min(self.s1.length, self.s2.length)
        for n in range(self.model.min_nucleation, num):
            for i, e1 in enumerate(self.nwise(self.s1.sequence, n)):
                for j, e2 in enumerate(self.nwise(self.s2.sequence, n)):
                    e1 = self.u(*e1)
                    e2 = self.u(*e2)
                    if self.model.space_dimensionality == '2D':
                        if i+j == len(ss)-n:
                            state = 'on_nucleation' if n == self.model.min_nucleation else 'zipping'
                        else: continue
                    else: 
                        if i+j != len(ss)-n: # --> Condition for checking that the nucleation is off register 
                            state = 'off_nucleation' if n == self.model.min_nucleation else 'backfray' 
                        else: state = 'on_nucleation' if n == self.model.min_nucleation else 'zipping'
                    if e1 == self.wc(e2[::-1]): 
                        if verbose: print(self.sab(self.s1.sequence, self.s2.sequence))
                        spacei = ' '*(i)
                        spacej = ' '*(j - i + len(ss) - n + 1)
                        if verbose: print(spacei+spacej.join([e1,e2]))
                        dpxdist = len(ss) - n - i - j
                        l = self.addpar(ss, i, n, '(')
                        r = self.addpar(ss, j, n, ')')
                        trap = self.sab(l,r)
                        obj = Complex(self.model, self.s1, self.s2, state=state, structure = trap, dpxdist=dpxdist)
                        self.DG.add_node(   trap,
                                            obj = obj, 
                                            pairs = obj.total_nucleations,
                                            state = obj.state,
                                            dpxdist = obj.dpxdist,
                                            tdx =(i,j),
                                            fre = obj.structureG())
                        dgtrap = obj.G
                        if n == self.model.min_nucleation:
                            if self.model.space_dimensionality == '3D':
                                fwd, bwd = self.nmethod(dgtrap)
                            elif self.model.space_dimensionality == '2D':
                                print('\n')
                                print(obj.structure)
                                print(obj.sdist)
                                fwd, bwd = self.nmethod(dgtrap, obj.sdist) #sdist is a measure of how many base pairs the current nucleation is away from the surface
                            self.DG.add_edge(self.simplex, trap, k = fwd, state = state)
                            self.DG.add_edge(trap, self.simplex, k = bwd, state = state)
                        elif n > self.model.min_nucleation:
                            self.sldbranches.append(dpxdist)
                            filter = self.filternodes('tdx', lambda x: (i-1 <= x[0] <= i+1) and (j-1 <= x[1] <= j+1), 
                                     self.filternodes('pairs', lambda x: x == n-1, 
                                     self.filternodes('dpxdist', lambda x: x == obj.dpxdist, self.DG)))
                            for node in filter.nodes():
                                if verbose: print(node, trap)
                                dgnode = self.DG.nodes[node]['fre']
                                fwd, bwd = self.zmethod('zipping', dgnode, dgtrap)
                                self.DG.add_edge(node, trap, k = fwd, state = 'backfray')
                                self.DG.add_edge(trap, node, k = bwd, state = 'backfray') 
                            if n == num-1: 
                                dgduplex = self.DG.nodes[self.duplex]['fre']
                                fwd, bwd = self.zmethod('zipping', dgtrap, dgduplex)
                                self.DG.add_edge(trap, self.duplex, k = fwd, state = 'zipping')
                                self.DG.add_edge(self.duplex, trap, k = bwd, state = 'zipping') 
        #get unique ids for all branches except branch 0 corresponding to on register nucleation
        self.sldbranches = set(self.sldbranches)

##########################################################################
##########################################################################

    def connect_slidings(self, verbose=False):
        for brc in combinations(self.sldbranches,2):
            # TODO add here a routine to insert missing slidings (due to lack of logic in get graph function)
            # connect slidings with eachother 
            if brc[0] != 0:
                leaf1 = self.filternodes('dpxdist', lambda x: x == brc[0], self.DG)
                mostables1 = []
                for comp in nx.strongly_connected_components(leaf1):
                    if len(comp) > 1:
                        subleaf = nx.subgraph(self.DG, list(comp))
                        mostables1.append(Structure(list(self.filternodes('fre', min, subleaf))[0]))
                leaf2 = self.filternodes('dpxdist', lambda x: x == brc[1], self.DG)
                mostables2 = []
                for comp in nx.strongly_connected_components(leaf2):
                    if len(comp) > 1:
                        subleaf = nx.subgraph(self.DG, list(comp))
                        mostables2.append(Structure(list(self.filternodes('fre', min, subleaf))[0]))
                for most1 in mostables1:
                    for most2 in mostables2:
                        self.DG.nodes[most1.str]['state'] = 'sliding'
                        self.DG.nodes[most2.str]['state'] = 'sliding'
                        if np.sign(brc[0])!=np.sign(brc[1]):
                            orig1, dest1, pkcond1 = self.kinetics.pkcond(most1, most2)
                            orig2, dest2, pkcond2 = self.kinetics.pkcond(most2, most1)
                            if pkcond1:
                                # pseudoknotting routine
                                #TODO add rates
                                if verbose: print(orig1.str,orig1.pktail_l,orig1.pktail_r,'-->',dest1.str,dest1.pktail_l,dest1.pktail_r,'pseudoknotting',wormdist)
                                self.DG.add_edge(orig1.str, dest1.str, k = 0, state = 'pseudoknotting')
                            if pkcond2:
                                if orig1 == orig2: pass
                                else:
                                    if verbose: print(orig2.str,orig2.pktail_l,orig2.pktail_r,'-->',dest2.str,dest2.pktail_l,dest2.pktail_r,'pseudoknotting back',wormdist)
                                    self.DG.add_edge(orig2.str, dest2.str, k = 0, state = 'pseudoknotting')                  
                            else: pass
                        else: 
                            # inchworming between slidings routine
                            wormdist = abs(brc[0]-brc[1])
                            #TODO add rates
                            # if verbose: print(most1.str,'-->',most2.str, 'inchworming', wormdist)
                            self.DG.add_edge(most1.str, most2.str, k = 0, state = 'inchworming')
                            self.DG.add_edge(most2.str, most1.str, k = 0, state = 'inchworming') 
            else:
                # inchworming towards duplex
                wormdist = abs(brc[1])
                most2 = Structure(list(
                    self.filternodes('fre', min,
                    self.filternodes('dpxdist', lambda x: x == brc[1], self.DG)
                    ))[0])
                #TODO add rates
                # print(most2.str,'-->',self.duplex,'inchworming duplex', wormdist)
                self.DG.add_edge(most2.str, self.duplex, k = 0, state = 'inchworming')  


##########################################################################

    def filternodes(self, property, function, graph):
        def filternode(node):
            try: 
                try:
                    return function(graph.nodes[node][property])
                except KeyError: pass
            except TypeError: 
                try: 
                    value = function([e[1][property] for e in list(graph.nodes.data())])
                    return graph.nodes[node][property] == value
                except KeyError: pass
        Rgraph = nx.subgraph_view(graph, filter_node=filternode)
        return Rgraph
    
    def filteredges(self, property, function, graph):
        def filteredge(node1, node2):
            try: 
                try:
                    return function(graph[node1][node2][property])
                except KeyError: pass
            except TypeError: 
                try: 
                    value = function([e[1][property] for e in list(graph.nodes.data())])
                    return graph[node1][node2][property] == value
                except KeyError: pass
        Rgraph = nx.subgraph_view(graph, filter_edge=filteredge)
        return Rgraph
    
##########################################################################

    def save_graph(self, PATH):
        import os 
        #convert node object to string of object type
        DGsave = self.DG.copy()
        for n in DGsave.nodes.data():
            n[1]['obj'] = str(type(n[1]['obj']))
            n[1]['fre'] = f"{n[1]['fre']:.3f}"
            try: n[1]['tdx'] = str(type(n[1]['tdx']))
            except KeyError: pass
        try: os.makedirs(PATH)
        except FileExistsError: pass 
        nx.write_gexf(DGsave,f'{PATH}/{self.s1.sequence}_graph_K.gexf')

###################
############## PROPERTIES
###################
    @property
    def simplex(self):
        return '+'.join(['.'*self.s1.length, '.'*self.s2.length])
    @property
    def duplex(self):
        return '+'.join(['('*self.s1.length, ')'*self.s2.length])
    @property
    def nodes(self):
        return self.DG.nodes()
    @property
    def displaysab(self):
        return '+'.join([self.s1.sequence,self.s2.sequence])

###################
############## UTILITIES
###################
########################################################
    def addpar(self, string, i, n, char):
        return string[:i] + char*n + string[i+n:]
########################################################
    def wc(self,a):
        wc = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join([wc[i] for i in a])
########################################################
    def ss(self,a):
        return '.'*len(a)
########################################################
    def sab(self,*args):
        return '+'.join(args)
########################################################
    def u(self,*args):
        return ''.join(args)
########################################################
    def nwise(self,iterable,n):
        iterators = tee(iterable, n)
        for i, iter in enumerate(iterators):
            for _ in range(i):
                next(iter, None)
        return zip(*iterators)

class FatalError(Exception):
    def __init__(self, message):
        super().__init__(message)


       # connect slidings with duplex 
        # duplex = Structure(self.duplex)   
        # for branch in self.sldbranches:
        #     mostable = list(
        #         self.filternodes('fre', min,
        #         self.filternodes('dpxdist', lambda x: x == branch, self.DG)
        #         ))[0]
        #     self.DG.nodes[mostable]['state'] = 'sliding'
        #     dgsliding = self.DG.nodes[mostable]['fre']
        #     # dgduplex = self.DG.nodes[self.duplex]['fre']
        #     fwd, _ = self.smethod('sliding', 0, dgsliding)                    
        #     if verbose: 
        #         dgstring = '{:.3f}'.format(dgsliding)
        #         fwdformat = '{:.3e}'.format(fwd)
        #         bwdformat = '{:.3e}'.format(0)
        #         # print(mostable, dgstring, self.kinetics.gammasliding(dgsliding))
        #         print(mostable, fwdformat, bwdformat, dgstring)
        #     #fwd = fwd / self.kinetics.gammasliding(dgsliding)# / abs(np.power(branch,1)) 
        #     #bwd = bwd / self.kinetics.gammasliding(dgsliding)# / abs(np.power(branch,1))
        #     self.DG.add_edge(mostable, self.duplex, k = fwd, state = 'sliding')
        #     self.DG.add_edge(self.duplex, mostable, k = 0, state = 'sliding')
        # for branch in self.sldbranches:
        #     leaf = self.filternodes('dpxdist', lambda x: x == branch, self.DG)
        #     components = nx.connected_components(leaf.to_undirected())
        #     for component in components:
        #         if len(component) > 1:
        #             subleaf = nx.subgraph(self.DG, list(component))
        #             mostable = list(self.filternodes('fre', min, subleaf).nodes())[0]
        #             self.DG.nodes[mostable]['state'] = 'sliding'
        #             dgsliding = self.DG.nodes[mostable]['fre']
        #             # dgduplex = self.DG.nodes[self.duplex]['fre']
        #             fwd, _ = self.smethod('sliding', 0, dgsliding)                    
        #             if verbose: 
        #                 dgstring = '{:.3f}'.format(dgsliding)
        #                 fwdformat = '{:.3e}'.format(fwd)
        #                 bwdformat = '{:.3e}'.format(0)
        #                 # print(mostable, dgstring, self.kinetics.gammasliding(dgsliding))
        #                 print(mostable, fwdformat, bwdformat, dgstring)
        #             #fwd = fwd / self.kinetics.gammasliding(dgsliding)# / abs(np.power(branch,1)) 
        #             #bwd = bwd / self.kinetics.gammasliding(dgsliding)# / abs(np.power(branch,1))
        #             self.DG.add_edge(mostable, self.duplex, k = fwd, state = 'sliding')
        #             self.DG.add_edge(self.duplex, mostable, k = 0, state = 'sliding')





                        # most1 = Structure(list(
                #     self.filternodes('fre', min,
                #     self.filternodes('dpxdist', lambda x: x == brc[0], self.DG)
                #     ))[0])
                # most2 = Structure(list(
                #     self.filternodes('fre', min,
                #     self.filternodes('dpxdist', lambda x: x == brc[1], self.DG)
                #     ))[0])