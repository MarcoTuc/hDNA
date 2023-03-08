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
                pass
                # # pick slidings in first branch
                # leaf1 = self.filternodes('dpxdist', lambda x: x == brc[0], self.DG)
                # mostables1 = []
                # for comp in nx.strongly_connected_components(leaf1):
                #     if len(comp) > 1:
                #         subleaf = nx.subgraph(self.DG, list(comp))
                #         mostables1.append(Structure(list(self.filternodes('fre', min, subleaf))[0]))
                # # pick slidings in second branch 
                # leaf2 = self.filternodes('dpxdist', lambda x: x == brc[1], self.DG)
                # mostables2 = []
                # for comp in nx.strongly_connected_components(leaf2):
                #     if len(comp) > 1:
                #         subleaf = nx.subgraph(self.DG, list(comp))
                #         mostables2.append(Structure(list(self.filternodes('fre', min, subleaf))[0]))
                # # iterate over all combinations 
                # for most1 in mostables1:
                #     for most2 in mostables2:
                #         # assign sliding state
                #         self.DG.nodes[most1.str]['state'] = 'sliding'
                #         self.DG.nodes[most2.str]['state'] = 'sliding'
                #         # check pseudoknotting conditions
                #         orig1, dest1, pkcond1 = self.kinetics.pkcond(most1, most2)
                #         orig2, dest2, pkcond2 = self.kinetics.pkcond(most2, most1)
                #         # assume pseudoknotting dominates over inchworming. 
                #         # e.g. if pseudo -> pseudo and only if not -> inchworm 
                        
                #         # TODO add a two-base pair formation routine like I did with the inchworming
                #         # and also weight the pseudoknotting transition tendency by using the relative
                #         # free energy difference since there is base pairing competition going on between
                #         # different registers in a pseudoknot 
                        
                #         if pkcond1:
                #             # pseudoknotting routine    
                #             self.psedoknotting_transition(orig1, dest1)                 
                #             # pkrate = self.kinetics.pkrate(orig1, dest1)
                #             if verbose: 
                #                 print(orig1.str,'-->',dest1.str,'pseudoknotting')
                #                 # print(f'{pkrate:.3e}')
                #             # self.DG.add_edge(orig1.str, dest1.str, k = pkrate, state = 'pseudoknotting')
                #         if pkcond2:
                #             if pkcond1:
                #                 if orig1.str == orig2.str: continue # don't duplicate 
                #             else:
                #                 self.psedoknotting_transition(orig2, dest2)
                #                 # pkrate = self.kinetics.pkrate(orig2, dest2)
                #                 if verbose: 
                #                     print(orig2.str,'-->',dest2.str,'pseudoknotting back')
                #                     # print(f'{pkrate:.3e}')
                #                 # self.DG.add_edge(orig2.str, dest2.str, k = pkrate, state = 'pseudoknotting')                  
                #         if (not pkcond1) and (not pkcond2): 
                #             # inchworming between slidings routine
                #             if most2.totbp > most1.totbp:
                #                 if most1.totbp > 2:
                #                     self.bulge_transition(most1, most2)
                #                     if verbose: print(most1.str,'-->',most2.str, 'inchworming')       
            # Connect slidings with duplex
            else:
                # pick slidings in second branch 
                leaf2 = self.filternodes('dpxdist', lambda x: x == brc[1], self.DG)
                mostables2 = [] 
                for comp in nx.strongly_connected_components(leaf2):
                    if len(comp) > 1:
                        subleaf = nx.subgraph(self.DG, list(comp))
                        mostables2.append(Structure(list(self.filternodes('fre', min, subleaf))[0]))
                # Iterate over most stable nodes for each leaf of the branch 
                for most2 in mostables2:
                    self.DG.nodes[most2.str]['state'] = 'sliding'
                    # pkcond = self.kinetics.pkconduplex(most2)
                    prob, pkcond = self.kinetics.pkcondsphereduplex(most2)
                    if pkcond: # pseudoknotting towards duplex 
                        self.psedoknotting_transition(most2, Structure(self.duplex), prob)
                        # pkrate = self.kinetics.pkrate(most2, Structure(self.duplex))
                        if verbose: 
                            print(most2.str,'-->',self.duplex,'pk duplex')
                            # print(f'{pkrate:.3e}')
                        # self.DG.add_edge(most2.str, self.duplex, k = pkrate, state = 'pseudoknotting')
                    else: # inchworming towards duplex
                        self.bulge_transition(most2, Structure(self.duplex))
    
    #OLD ROUTINE
    def connect_slidings(self, verbose=True):
        self.sldbranches.remove(0)
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
                    gsl    = self.kinetics.gammasliding(dgsliding)       
                    if verbose: 
                        dgstring = '{:.3f}'.format(dgsliding)
                        fwdformat = '{:.3e}'.format(fwd)
                        bwdformat = '{:.3e}'.format(0)
                        # print(mostable, dgstring, self.kinetics.gammasliding(dgsliding))
                        print(mostable, fwdformat, bwdformat, dgstring)
                    #fwd = fwd / self.kinetics.gammasliding(dgsliding)# / abs(np.power(branch,1)) 
                    #bwd = bwd / self.kinetics.gammasliding(dgsliding)# / abs(np.power(branch,1))
                    self.DG.add_edge(mostable, self.duplex, k = fwd, state = 'sliding')
                    self.DG.add_edge(self.duplex, mostable, k = 0, state = 'sliding')


    def bulge_transition(self, origin, target):
        willbulge = list(self.DG.neighbors(origin.str))
        try: willbulge.remove(target)
        except ValueError: pass 

        for e in willbulge[::-1]:
            # print(e[0])
            if e[0] == 'B':
                # print('yo', e)
                willbulge.remove(e)
            elif e[0] == 'P':
                # print('yo', e)
                willbulge.remove(e)
            else: continue
            
        neighs = [Structure(e) for e in willbulge]
        # inchworming happens from the neighbors of the fully slided state
        # and not from the slided state itself, so we connect these to the target
        for wb in neighs:
            # create an abstract bulged state
            bulge = f'BULGED_{wb.str}'
            self.DG.add_node(bulge, 
                            obj = '',
                            pairs = origin.totbp, 
                            state = 'bulgeloop',
                            dpxdist = f'{origin.register} - {target.register}',
                            fre = '')
            kbulged = self.kinetics.bulging(origin, target)
            self.DG.add_edge(wb.str, bulge, k = kbulged, state='bulging')
            # self.DG.add_edge(bulge, wb.str, k = self.kinetics.avgunzip(), state='bulging')
            # create a second bulging state, if reached proceed with irreversible bulging transition
            # The second bulging happens at speed equal to the first bulging kinda 
            bulging = f'BULGING_{wb.str}'
            self.DG.add_node(bulging, 
                            obj = '',
                            pairs = origin.totbp, 
                            state = 'bulgeloop',
                            dpxdist = f'{origin.register} - {target.register}',
                            fre = '')
            kbulging = kbulged #self.kinetics.bulging(origin, target)
            # tau = (1/kbulging) + (1/self.kinetics.avgunzip())
            self.DG.add_edge(bulge, bulging, k = kbulging, state = 'bulging')
            kinchworm = kbulged #self.kinetics.iwrate(origin, target)
            print('inchworm')
            print(origin.str, target.str)
            print('kbulged ', f'{kbulged:.3e}')
            print('kbulging',f'{kbulging:.3e}')
            print('inchworm',f'{kinchworm:.3e}')
            print()
            self.DG.add_edge(bulging, target.str, k = kinchworm, state='inchworming')

    def psedoknotting_transition(self, origin, target, probability): 
        # combinatorics = 1/(2*(origin.tail_ll + origin.tail_lr))
        print('pseudoknot probability:',probability)
        # print(combinatorics)
        pseudoknot = f'PKCOLL_{origin.str}'
        self.DG.add_node(pseudoknot, 
                         obj = '',
                         pairs = origin.totbp, 
                         state = 'pseudoknot',
                         dpxdist = f'{origin.register} - {target.register}',
                         fre = '')
        pkcollrate = self.kinetics.nucleationrate * probability
        self.DG.add_edge(origin.str, pseudoknot, k=pkcollrate, state='pseudoknotting')
        freorigin = self.DG.nodes[origin.str]['fre']
        fretarget = self.DG.nodes[target.str]['fre']
        localZ = self.kinetics.localZ(freorigin, fretarget-freorigin)
        Zorigin = self.kinetics.genboltz(freorigin)/localZ
        Ztarget = self.kinetics.genboltz(fretarget-freorigin)/localZ
        pktransratefwd = self.kinetics.pkrate(origin, target) * abs(Ztarget)
        pktransratebwd = self.kinetics.pkrate(target, origin) * abs(Zorigin)
        # pktransratefwd = pktransratefwd*self.kinetics.genboltz(freorigin)
        print('pseudoknot')
        print(origin.str, target.str)
        print('collision', f'{pkcollrate:.3e}')
        print('fwd:',f'{pktransratefwd:.3e}','bwd:', f'{pktransratebwd:.3e}')
        print('fres',freorigin, fretarget)
        print('zeta',Ztarget, Zorigin)
        print()
        # fwd, _ = self.kinetics.kawasaki('any', freorigin, fretarget, rate=pktransratefwd) * frebalance
        # bwd, _ = self.kinetics.kawasaki('any', fretarget, freorigin, rate=pktransratebwd) / frebalance
        self.DG.add_edge(pseudoknot, target.str, k = pktransratefwd, state = 'pseudoknotting')
        self.DG.add_edge(pseudoknot, origin.str, k = pktransratebwd, state = 'pseudoknotting')



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
            try: n[1]['fre'] = f"{n[1]['fre']:.3f}"
            except ValueError: pass
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



                # willbulge = list(self.DG.neighbors(most2.str))
                        # try: willbulge.remove(self.duplex)
                        # except ValueError: pass 
                        # neighs = [Structure(e) for e in willbulge]
                        # # inchworming happens from the neighbors of the fully slided state
                        # # and not from the slided state itself, so we connect these to the duplex
                        # for wb in neighs:
                        #     # create an abstract bulged state
                        #     bulge = f'BULGED_{wb.str}'
                        #     self.DG.add_node(bulge, 
                        #                     obj = '',
                        #                     pairs = most2.totbp, 
                        #                     state = 'bulgeloop',
                        #                     dpxdist = f'{most2.register} - {0}',
                        #                     fre = '')
                            
                        #     kbulging = self.kinetics.bulging(most2, Structure(self.duplex))
                        #     self.DG.add_edge(wb.str, bulge, k = kbulging, state='bulging')
                        #     self.DG.add_edge(bulge, wb.str, k = self.kinetics.avgunzip(), state='bulging')
                            
                        #     # create a second bulging state, if reached proceed with irreversible bulging transition
                        #     # The second bulging happens at speed equal to the first bulging kinda 
                        #     bulging = f'BULGING_{most2.str}'
                        #     self.DG.add_node(bulging, 
                        #                     obj = '',
                        #                     pairs = most2.totbp, 
                        #                     state = 'bulgeloop',
                        #                     dpxdist = most2.register,
                        #                     fre = '')
                        #     kbulging = self.kinetics.bulging(most2, Structure(self.duplex))
                        #     tau = (1/kbulging) + (1/self.kinetics.avgunzip())
                        #     self.DG.add_edge(bulge, bulging, k = 1/tau, state = 'bulging')

                        #     kinchworm = self.kinetics.iwrate(most2)
                        #     self.DG.add_edge(bulging, self.duplex, k = kinchworm, state='inchworming')
                        
                        
                        # most2 = Structure(list(
                        #     self.filternodes('fre', min,
                        #     self.filternodes('dpxdist', lambda x: x == brc[1], self.DG)
                        #     ))[0])