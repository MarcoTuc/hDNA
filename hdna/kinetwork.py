import networkx as nx
import numpy as np 
import pandas as pd 
import nupack as nu 
from itertools import pairwise, combinations, tee

from .complex import Complex
from .model import Model
from .strand import Strand
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
        print(type(self.dxobj.G))
        self.f2exit = []
        self.f1exit = []
        self.get_graphmid()
        # self.connect_slidings()


##########################################################################
##########################################################################
    
    def get_graphleven(self, verbose=False):

        self.sldbranches = []
        ss = self.simplex.split('+')[0] 
        minlen = min(self.s1.length, self.s2.length)
        
        for n in range(self.model.min_nucleation, minlen):
            for i, e1 in enumerate(self.nwise(self.s1.sequence, n)):
                for j, e2 in enumerate(self.nwise(self.s2.sequence, n)):
                    
                    e1 = self.u(*e1)
                    e2 = self.u(*e2)
                    
                    # if self.model.space_dimensionality == '2D':
                    #     if i+j == len(ss)-n:
                    #         state = 'on_nucleation' if n == self.model.min_nucleation else 'zipping'
                    #     else: continue
                    # else: 
                    #     if i+j != len(ss)-n: # --> Condition for checking that the nucleation is off register 
                    #         state = 'off_nucleation' if n == self.model.min_nucleation else 'backfray' 
                    #     else: state = 'on_nucleation' if n == self.model.min_nucleation else 'zipping'
                    
                    newbps = self.wcwhere(e1, e2[::-1])
                    print(newbps)
### IDEA
# Make all the one base pair long nucleations
# All possible subsequent backfrays will be compositions of nucleations at the same 
# duplexdistance/register. Then use levenshtein distance between structures to determine
# For each nucleation in each register what's the next guy coming up. Connect these. 
# This approach should work and should also be flexible. 


                    if sum(newbps)==n: 

                        if verbose: 
                            print(self.sab(self.s1.sequence, self.s2.sequence))
                            spacei = ' '*(i)
                            spacej = ' '*(j - i + len(ss) - n + 1)
                            print(spacei+spacej.join([e1,e2]))

                        dpxdist = len(ss) - n - i - j

                        l = self.addparwhere(ss, i, newbps, '(')
                        r = self.addparwhere(ss, j, newbps[::-1], ')')

                        trap = self.sab(l,r)
                        print(self.sab(self.s1.sequence, self.s2.sequence))
                        print(trap, dpxdist)
        #                 obj = Complex(self.model, self.s1, self.s2, state=state, structure = trap, dpxdist=dpxdist)
        #                 self.DG.add_node(   trap,
        #                                     obj = obj, 
        #                                     pairs = obj.total_nucleations,
        #                                     state = obj.state,
        #                                     dpxdist = obj.dpxdist,
        #                                     tdx =(i,j),
        #                                     fre = obj.structureG())
        #                 dgtrap = obj.G
        #                 if n == self.model.min_nucleation:
        #                     if self.model.space_dimensionality == '3D':
        #                         fwd, bwd = self.nmethod(dgtrap)
        #                     elif self.model.space_dimensionality == '2D':
        #                         if verbose: 
        #                             print('\n')
        #                             print(obj.structure)
        #                             print(obj.sdist)
        #                         fwd, bwd = self.nmethod(dgtrap, obj.sdist) #sdist is a measure of how many base pairs the current nucleation is away from the surface
        #                     self.DG.add_edge(self.simplex, trap, k = fwd, state = state)
        #                     self.DG.add_edge(trap, self.simplex, k = bwd, state = state)
        #                 elif n > self.model.min_nucleation:
        #                     self.sldbranches.append(dpxdist)
        #                     f1 = self.filternodes('dpxdist', lambda x: x == obj.dpxdist, self.DG)
        #                     f2 = self.filternodes('pairs', lambda x: x == n-1, f1)
        #                     f3 = self.filternodes('tdx', lambda x: (i-1 <= x[0] <= i+1) and (j-1 <= x[1] <= j+1), f2)
        #                     print(*list(f3.nodes()),sep='\n')
        #                     print(obj.dpxdist)
        #                     for node in f3.nodes():
        #                         if node != trap:
        #                             # print('f3has',node)
        #                             if verbose: print(node, trap)
        #                             dgnode = self.DG.nodes[node]['fre']
        #                             fwd, bwd = self.zmethod('zipping', dgnode, dgtrap)
        #                             # print(node, trap)
        #                             self.DG.add_edge(node, trap, k = fwd, state = 'backfray')
        #                             self.DG.add_edge(trap, node, k = bwd, state = 'backfray') 
        #                     if n == minlen-1:
        #                         print('yo', trap) 
        #                         dgduplex = self.DG.nodes[self.duplex]['fre']
        #                         fwd, bwd = self.zmethod('zipping', dgtrap, dgduplex)
        #                         self.DG.add_edge(trap, self.duplex, k = fwd, state = 'zipping')
        #                         self.DG.add_edge(self.duplex, trap, k = bwd, state = 'zipping') 
        # #get unique ids for all branches except branch 0 corresponding to on register nucleation
        # self.sldbranches = set(self.sldbranches)
        # self.sldbranches.remove(0)
    
    def get_graphnew(self, verbose=False):

        self.sldbranches = []
        ss = self.simplex.split('+')[0] 
        branches = []

        for i, e1 in enumerate(self.s1.sequence):
            for j, e2 in enumerate(self.s2.sequence):
                    e1 = self.u(*e1)
                    e2 = self.u(*e2)
                    newbps = self.wcwhere(e1, e2[::-1])
                    if any(newbps):
                        dpxdist = len(ss) - 1 - i - j
                        l = self.addparwhere(ss, i, newbps, '(')
                        r = self.addparwhere(ss, j, newbps[::-1], ')')
                        trap = self.sab(l,r)
                        print(self.sab(self.s1.sequence, self.s2.sequence))
                        print(trap, dpxdist)
                        branches.append({trap: dpxdist})
        return branches



    def get_graphmid(self, verbose=False):

        self.sldbranches = []
        ss = self.simplex.split('+')[0] 
        minlen = min(self.s1.length, self.s2.length)
        
        for n in range(self.model.min_nucleation, minlen):
            for i, e1 in enumerate(self.nwise(self.s1.sequence, n)):
                for j, e2 in enumerate(self.nwise(self.s2.sequence, n)):
                    
                    e1 = self.u(*e1)
                    e2 = self.u(*e2)
                    
                    # if self.model.space_dimensionality == '2D':
                    #     if i+j == len(ss)-n:
                    #         state = 'on_nucleation' if n == self.model.min_nucleation else 'zipping'
                    #     else: continue
                    # else: 
                    #     if i+j != len(ss)-n: # --> Condition for checking that the nucleation is off register 
                    #         state = 'off_nucleation' if n == self.model.min_nucleation else 'backfray' 
                    #     else: state = 'on_nucleation' if n == self.model.min_nucleation else 'zipping'
                    
                    newbps = self.wcwhere(e1, e2[::-1])
    
                    if sum(newbps)==n: 

                        if verbose: 
                            print(self.sab(self.s1.sequence, self.s2.sequence))
                            spacei = ' '*(i)
                            spacej = ' '*(j - i + len(ss) - n + 1)
                            print(spacei+spacej.join([e1,e2]))

                        dpxdist = len(ss) - n - i - j

                        l = self.addparwhere(ss, i, newbps, '(')
                        r = self.addparwhere(ss, j, newbps[::-1], ')')

                        trap = self.sab(l,r)
                        print(self.sab(self.s1.sequence, self.s2.sequence))
                        print(trap, dpxdist)
        #                 obj = Complex(self.model, self.s1, self.s2, state=state, structure = trap, dpxdist=dpxdist)
        #                 self.DG.add_node(   trap,
        #                                     obj = obj, 
        #                                     pairs = obj.total_nucleations,
        #                                     state = obj.state,
        #                                     dpxdist = obj.dpxdist,
        #                                     tdx =(i,j),
        #                                     fre = obj.structureG())
        #                 dgtrap = obj.G
        #                 if n == self.model.min_nucleation:
        #                     if self.model.space_dimensionality == '3D':
        #                         fwd, bwd = self.nmethod(dgtrap)
        #                     elif self.model.space_dimensionality == '2D':
        #                         if verbose: 
        #                             print('\n')
        #                             print(obj.structure)
        #                             print(obj.sdist)
        #                         fwd, bwd = self.nmethod(dgtrap, obj.sdist) #sdist is a measure of how many base pairs the current nucleation is away from the surface
        #                     self.DG.add_edge(self.simplex, trap, k = fwd, state = state)
        #                     self.DG.add_edge(trap, self.simplex, k = bwd, state = state)
        #                 elif n > self.model.min_nucleation:
        #                     self.sldbranches.append(dpxdist)
        #                     f1 = self.filternodes('dpxdist', lambda x: x == obj.dpxdist, self.DG)
        #                     f2 = self.filternodes('pairs', lambda x: x == n-1, f1)
        #                     f3 = self.filternodes('tdx', lambda x: (i-1 <= x[0] <= i+1) and (j-1 <= x[1] <= j+1), f2)
        #                     print(*list(f3.nodes()),sep='\n')
        #                     print(obj.dpxdist)
        #                     for node in f3.nodes():
        #                         if node != trap:
        #                             # print('f3has',node)
        #                             if verbose: print(node, trap)
        #                             dgnode = self.DG.nodes[node]['fre']
        #                             fwd, bwd = self.zmethod('zipping', dgnode, dgtrap)
        #                             # print(node, trap)
        #                             self.DG.add_edge(node, trap, k = fwd, state = 'backfray')
        #                             self.DG.add_edge(trap, node, k = bwd, state = 'backfray') 
        #                     if n == minlen-1:
        #                         print('yo', trap) 
        #                         dgduplex = self.DG.nodes[self.duplex]['fre']
        #                         fwd, bwd = self.zmethod('zipping', dgtrap, dgduplex)
        #                         self.DG.add_edge(trap, self.duplex, k = fwd, state = 'zipping')
        #                         self.DG.add_edge(self.duplex, trap, k = bwd, state = 'zipping') 
        # #get unique ids for all branches except branch 0 corresponding to on register nucleation
        # self.sldbranches = set(self.sldbranches)
        # self.sldbranches.remove(0)


    def get_graphold(self, verbose=False):
        self.sldbranches = []
        ss = self.simplex.split('+')[0] 
        minlen = min(self.s1.length, self.s2.length)
        
        for n in range(self.model.min_nucleation, minlen):
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
                        
                        if verbose: 
                            print(self.sab(self.s1.sequence, self.s2.sequence))
                            spacei = ' '*(i)
                            spacej = ' '*(j - i + len(ss) - n + 1)
                            print(spacei+spacej.join([e1,e2]))

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
                                if verbose: 
                                    print('\n')
                                    print(obj.structure)
                                    print(obj.sdist)
                                fwd, bwd = self.nmethod(dgtrap, obj.sdist) #sdist is a measure of how many base pairs the current nucleation is away from the surface
                            self.DG.add_edge(self.simplex, trap, k = fwd, state = state)
                            self.DG.add_edge(trap, self.simplex, k = bwd, state = state)
                        elif n > self.model.min_nucleation:
                            self.sldbranches.append(dpxdist)
                            f1 = self.filternodes('dpxdist', lambda x: x == obj.dpxdist, self.DG)
                            self.f1exit.append(f1)
                            f2 = self.filternodes('pairs', lambda x: x == n-1, f1)
                            self.f2exit.append(f2)
                            f3 = self.filternodes('tdx', lambda x: (i-1 <= x[0] <= i+1) and (j-1 <= x[1] <= j+1), f2)
                            for node in f3.nodes():
                                if verbose: print(node, trap)
                                dgnode = self.DG.nodes[node]['fre']
                                fwd, bwd = self.zmethod('zipping', dgnode, dgtrap)
                                self.DG.add_edge(node, trap, k = fwd, state = 'backfray')
                                self.DG.add_edge(trap, node, k = bwd, state = 'backfray') 
                            if n == minlen-1:
                                print('yo', trap) 
                                dgduplex = self.DG.nodes[self.duplex]['fre']
                                fwd, bwd = self.zmethod('zipping', dgtrap, dgduplex)
                                self.DG.add_edge(trap, self.duplex, k = fwd, state = 'zipping')
                                self.DG.add_edge(self.duplex, trap, k = bwd, state = 'zipping') 
        #get unique ids for all branches except branch 0 corresponding to on register nucleation
        self.sldbranches = set(self.sldbranches)
        self.sldbranches.remove(0)

##########################################################################
##########################################################################

    def connect_slidings(self, verbose=True):
        for branch in self.sldbranches:
            leaf = self.filternodes('dpxdist', lambda x: x == branch, self.DG)
            components = nx.connected_components(leaf.to_undirected())
            for component in components:
                if len(component) > 1:
                    subleaf = nx.subgraph(self.DG, list(component))
                    mostable = list(self.filternodes('fre', min, subleaf).nodes())[0]
                    self.DG.nodes[mostable]['state'] = 'sliding'
                    dgsliding = self.DG.nodes[mostable]['fre']
                    # dgduplex = self.DG.nodes[self.duplex]['fre']
                    fwd, _ = self.smethod('sliding', 0, dgsliding)                    
                    if verbose: 
                        dgstring = '{:.3f}'.format(dgsliding)
                        fwdformat = '{:.3e}'.format(fwd)
                        bwdformat = '{:.3e}'.format(0)
                        print(mostable, fwdformat, bwdformat, dgstring)
                    self.DG.add_edge(mostable, self.duplex, k = fwd, state = 'sliding')
                    self.DG.add_edge(self.duplex, mostable, k = 0, state = 'sliding')

##########################################################################

    def filternodes(self, property, function, graph='self'):
        if graph == 'self':
            graph = self.DG
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
    def addparwhere(self, string, i, where, char):
        for j, wc in enumerate(where):
            if wc:
                string = string[:i+j] + char + string[i+j+1:]
        return string
    def addpariter(self, string, i, where, char):
        for j, wc in enumerate(where):
            if wc:
                string = string[:i+j] + char + string[i+j+1:]
                yield string
########################################################
    def wc(self,a):
        wc = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join([wc[i] for i in a])
    def wcwhere(self,a,b):
        if len(a) != len(b):
            raise ValueError('Cannot confront slices of different length')
        wcd = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return [a[i] == wcd[b[i]] for i in range(len(a))]

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








# I tried to improve the graph creation flow by modifying the iteration logic. 
# The attempt was not successful since it would've taken lots of time to become functional, 
# time which I didn't have, although I believe I was onto something here. 

# def get_graph(self, verbose=False):
#         self.sldbranches = []
#         ss = self.simplex.split('+')[0] 
#         lenreg = min(self.s1.length, self.s2.length)
#         # iterate over all registers 
#         for r in range(self.model.min_nucleation - lenreg, lenreg - self.model.min_nucleation + 1):
#             # number of base pairs in the register
#             nbpr = lenreg - abs(r)
#             print(r)
#             # slice sequences according to the register 
#             if np.sign(r) == -1:
#                 space1 = 0
#                 space2 = 0
#                 slice1 = self.s1.sequence[:nbpr]
#                 slice2 = self.s2.sequence[:nbpr]
#             else:
#                 space1 = lenreg-nbpr
#                 space2 = 0
#                 slice1 = self.s1.sequence[lenreg-nbpr:]
#                 slice2 = self.s2.sequence[:nbpr]
#             for n in range(nbpr):
#                 for i, e1 in enumerate(self.nwise(slice1, n)):
#                     for j, e2 in enumerate(self.nwise(slice2, n)):
#                         e1 = self.u(*e1)
#                         e2 = self.u(*e2)
#                         if self.model.space_dimensionality == '2D':
#                             if r == 0:
#                                 nstate = 'on_nucleation' if n == self.model.min_nucleation else 'zipping'
#                             else: continue
#                         else: 
#                             if r != 0: # --> Condition for checking that the nucleation is off register 
#                                 nstate = 'off_nucleation' if n == self.model.min_nucleation else 'backfray' 
#                             else: nstate = 'on_nucleation' if n == self.model.min_nucleation else 'zipping'
#                         print(e1,i,e2,j)
#                         newbps = self.wcwhere(e1, e2) #returns a list of true values where wc pairing is found 
#                         print(newbps)
#                         if any(newbps):
#                             print('spaces',space1, space2)
#                             leftiter = self.addpariter(ss, i+space1, newbps, '(')
#                             righiter = self.addpariter(ss, j+space2, newbps, ')')
#                             for left, righ in zip(leftiter, righiter):
#                                 trap = self.sab(left,righ)
#                                 print(self.s1.sequence+'+'+self.s2.sequence)
#                                 print(trap)
#                                 obj = Complex(self.model, self.s1, self.s2, state=nstate, structure = trap, dpxdist=r)
#                                 self.DG.add_node(   trap,
#                                                     obj = obj, 
#                                                     pairs = obj.total_nucleations,
#                                                     state = obj.state,
#                                                     dpxdist = obj.dpxdist,
#                                                     tdx =(i,j),
#                                                     fre = obj.G)
#                                 dgtrap = obj.G
#                                 # Add nucleations wether on core or off core 
#                                 if n == self.model.min_nucleation:
#                                     if self.model.space_dimensionality == '3D':
#                                         fwd, bwd = self.nmethod(dgtrap)
#                                     elif self.model.space_dimensionality == '2D':                    
#                                         fwd, bwd = self.nmethod(dgtrap, obj.sdist) #sdist is a measure of how many base pairs the current nucleation is away from the surface
#                                     self.DG.add_edge(self.simplex, trap, k = fwd, state = nstate)
#                                     self.DG.add_edge(trap, self.simplex, k = bwd, state = nstate)
#                                 # Add subsequent backfrays or zippings
#                                 elif n > self.model.min_nucleation:
#                                     self.sldbranches.append(r)
#                                     f1 = self.filternodes('dpxdist', lambda x: x == obj.dpxdist, self.DG)
#                                     # self.f1exit.append(f1)
#                                     f2 = self.filternodes('pairs', lambda x: x == n-1, f1)
#                                     # self.f2exit.append(f2)
#                                     f3 = self.filternodes('tdx', lambda x: (i-1 <= x[0] <= i+1) and (j-1 <= x[1] <= j+1), f2)
#                                     for node in f3.nodes():
#                                         state = 'backfray' if r != 0 else 'zipping'
#                                         dgnode = self.DG.nodes[node]['fre']
#                                         fwd, bwd = self.zmethod('zipping', dgnode, dgtrap)
#                                         self.DG.add_edge(node, trap, k = fwd, state = state)
#                                         self.DG.add_edge(trap, node, k = bwd, state = state) 
#                                     # Only useful to connect last zippings with duplex
#                                     if n == nbpr-1:
#                                         if r == 0:
#                                             dgduplex = self.DG.nodes[self.duplex]['fre']
#                                             fwd, bwd = self.zmethod('zipping', dgtrap, dgduplex)
#                                             self.DG.add_edge(trap, self.duplex, k = fwd, state = 'zipping')
#                                             self.DG.add_edge(self.duplex, trap, k = bwd, state = 'zipping') 
#                                         else: pass

#         #get unique ids for all branches except branch 0 corresponding to on register nucleation
#         self.sldbranches = set(self.sldbranches)
#         self.sldbranches.remove(0)