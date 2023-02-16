import networkx as nx
import numpy as np 
import pandas as pd 

from itertools import pairwise, combinations, tee

from .complex import Zippo, Sliding, Complex
from .chamber import Chamber
from .model import Model
from .strand import Strand
 
class Kinetwork(object):

    """
    Given a chamber it will generate the corresponding
    Kinetic network to be later passed to the simulator.
    To do so it composes two graphs: 
        - Zipping Graph 
        - Sliding Graph
    The Zipping Graph is the composition of single zipping graphs ---> See zippingraph() method
    corresponding to each oncore taken from the chamber.
    The Sliding Graph is the composition of single sliding graphs ---> See slidingraph() method
    """

    def __init__(self, model: Model, s1: Strand, s2: Strand, clean=False):

        self.model = model 
        self.s1 = s1        #53
        self.s2 = s2.invert        #35
        self.min_nucleation = self.model.min_nucleation
        self.sliding_cutoff = self.model.sliding_cutoff
        self.chamber = Chamber(self.model, self.s1, self.s2.invert)
        self.kinetics = Kinetics(self.model, self.chamber)
        self.kinetics.set_slidingrate(self.model.sliding)
        self.kinetics.set_zippingrate(self.model.zipping)
        self.simplex = self.chamber.singlestranded.structure
        self.duplex = self.chamber.duplex.structure

        self.nmethod = self.kinetics.kawasaki
        self.zmethod = self.kinetics.kawasaki
        self.smethod = self.kinetics.metropolis

        if not clean: 
            self.zippingraph()
            self.slidingraph()
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


    def completegraph(self):
        self.DG = nx.compose(self.SG, self.ZG)
        for node in self.DG.nodes():
            self.DG.nodes[node]['fre'] = self.DG.nodes[node]['obj'].G
        

    def zippingraph(self):
        self.ZG = nx.DiGraph()
        self.ZG.add_node(self.simplex, 
                         obj = self.chamber.singlestranded,
                         pairs = 0, 
                         state = 'singlestranded',
                         dpxdist = self.chamber.singlestranded.dpxdist)
        self.ZG.add_node(self.duplex, 
                         obj = self.chamber.duplex, 
                         pairs = self.chamber.duplex.total_nucleations, 
                         state = 'duplex',
                         dpxdist = 0)
        #general subfunction to compute zippingraphs of each oncore 
        def subzippingraph(oncore):
            #utility to get indices of string update
            def get_i(lst):
                indices = []
                for i, (el1, el2) in enumerate(pairwise(lst)):
                    if el1 == '.' and el2 != '.':
                        indices.append(i)
                    if el1 != '.' and el2 == '.':
                        indices.append(i+1)
                return indices
            #updates current the given structure
            def update_structure(string, side, character: str):
                updated = string
                moveslist = []
                moves = get_i(string)[::side]
                for m in moves:
                    up = updated[:m] + character + updated[m+1:]
                    moveslist.append(up)
                return moveslist
            #recursively add nodes and edges to their leafs in a tree like manner until duplexation
            #duplexation is automatically recognised by the production of an empty moveslist list 
            def leafs(left, right, graph):
                lmoves = update_structure(left, 1, '(')
                rmoves = update_structure(right, -1, ')')
                inbound = '+'.join([left, right])
                zobj = Zippo(self.model, self.s1, self.s2, state='zipping', structure=inbound, clean = True)
                graph.add_node( inbound,
                                obj = zobj,
                                state = zobj.state,
                                pairs = zobj.total_nucleations)
                for lmove, rmove in zip(lmoves, rmoves):
                    outbound = '+'.join([lmove,rmove])
                    zobj = Zippo(self.model, self.s1, self.s2, state='zipping', structure=outbound, clean = True)
                    graph.add_node(outbound,
                                   obj = zobj,
                                   state = zobj.state,
                                   pairs = zobj.total_nucleations)
                    graph.add_edge(inbound, outbound)
                    graph.add_edge(outbound, inbound)
                    yield list(leafs(lmove,rmove,graph))
            #takes inputted oncore and splits it to update structure simultaneously on the left and the right
            left, right = oncore.structure.split('+')
            subZG = nx.DiGraph()
            list(leafs(left, right, subZG))
            return subZG 
        #add the subzippingraph to the general zipping graph for each oncore
        for on in self.chamber.oncores:
            addgraph = subzippingraph(on)
            self.ZG.update(addgraph)

        for node in self.ZG.nodes():
            self.ZG.nodes[node]['obj'].structureG()
        
        nucleations = self.filternodes('pairs', lambda x: x == self.min_nucleation, self.ZG)
        for node in list(nucleations.nodes()):
            state = 'on_nucleation'
            self.ZG.nodes[node]['state'] = state
            dgss = 0 #reference simplex state
            dgnuc = self.ZG.nodes[node]['obj'].G #nucleation free energy
            # compute forward and backward rates
            fwd, bwd = self.nmethod(state, dgss, dgnuc) 
            nucnorm = (self.s1.length + self.s2.length - 1) * self.chamber.duplex.total_nucleations
            fwd = fwd/nucnorm
            bwd = bwd/nucnorm
            self.ZG.add_edge(self.simplex, node, k = fwd, state = state)
            self.ZG.add_edge(node, self.simplex, k = bwd, state = state)
            next = list(self.ZG.neighbors(node))
            next.remove(self.simplex)
            for e in next:
                dge = self.ZG.nodes[e]['obj'].G
                # compute forward and backward rates
                fwd, bwd = self.zmethod('zipping', dgnuc, dge) 
                self.ZG.add_edge(node, e, k = fwd, state = 'zipping')
                self.ZG.add_edge(e, node, k = bwd, state = 'zipping')
                
        
        mbare = list(set(self.ZG.nodes()) - set(nucleations.nodes()))
        for e1, e2 in nx.subgraph_view(self.ZG, lambda x: x in mbare).to_undirected().edges(): 
            dg1 = self.ZG.nodes[e1]['obj'].G
            dg2 = self.ZG.nodes[e2]['obj'].G
            fwd, bwd = self.zmethod('zipping', dg1, dg2)
            self.ZG[e1][e2]['k'] = fwd
            self.ZG[e2][e1]['k'] = bwd
            self.ZG[e1][e2]['state'] = 'zipping' 
            self.ZG[e2][e1]['state'] = 'zipping'

        self.ZG.nodes[self.duplex]['state'] = 'duplex'

##############################################################################################            
    
    def slidingraph(self):
        self.SG = nx.DiGraph()

        self.SG.add_node(self.simplex, 
                         obj = self.chamber.singlestranded,
                         pairs = 0, 
                         state = 'singlestranded',
                         dpxdist = self.chamber.singlestranded.dpxdist)
        self.SG.add_node(self.duplex, 
                         obj = self.chamber.duplex, 
                         state = 'duplex',
                         pairs = self.chamber.duplex.total_nucleations, 
                         dpxdist = 0)

        def backfraygraph(sliding, side):
            #get indices where the sliding has possible basepairs by parsing its completed structure
            def getix(string):
                lix = []
                rix = []
                for i, e in enumerate(list(string)):
                    if e == '(':
                        lix.append(i)
                    if e == ')':
                        rix.append(i)
                return [(l,r) for l,r in zip(lix, rix[::-1])]
            #confront currently given structure with target sliding to get adjacent next basepairs to form
            def adjacent(source, target):
                adj = []
                tx = getix(target)
                sx = getix(source)
                for i, item in enumerate(tx):
                    if len(sx) > 1:
                        check = [sx[0], sx[-1]]
                        if item in check:
                            if i == 0:
                                pass
                                # adj.append(tx[i+1])
                            elif i == len(tx)-1:
                                pass
                                # adj.append(tx[i-1])
                            else:
                                if tx[i-1] not in sx:
                                    if tx[i-1] not in adj:
                                        adj.append(tx[i-1])
                                if tx[i+1] not in sx:
                                    if tx[i+1] not in adj:
                                        adj.append(tx[i+1])
                    else:
                        check = sx
                        if item in check:
                            if i == 0:

                                adj.append(tx[i+1])
                            elif i == len(tx)-1:
                                adj.append(tx[i-1])
                            else:
                                adj.append(tx[i-1])
                                adj.append(tx[i+1])
                return adj 
            #update once every new adjacent basepair to form
            def update_structure(structure, target):
                updated = structure
                moveslist = []
                moves = adjacent(structure, target)
                for m in moves:
                    up = updated[:m[0]] + '(' + updated[m[0]+1:m[1]] + ')' + updated[m[1]+1:]
                    moveslist.append(up)
                return moveslist
            #leafing tree recursive function for generating all moves 
            def leafs(structure, target, graph):
                moves = update_structure(structure, target)
                obj = Complex(  self.model,
                                self.s1,
                                self.s2, 
                                state='backfray',
                                structure = structure,
                                dpxdist=sliding.dpxdist,
                                clean = True)
                obj.dtc = sliding.total_nucleations - obj.total_nucleations
                obj.side = side
                graph.add_node( structure,
                                obj = obj, 
                                pairs = obj.total_nucleations,
                                state = obj.state,
                                dpxdist = obj.dpxdist,
                                dtc = obj.dtc,
                                side = obj.side)
                for move in moves:
                    obj = Complex(
                                self.model,
                                self.s1,
                                self.s2, 
                                state='backfray',
                                structure = move,
                                dpxdist=sliding.dpxdist,
                                clean = True)
                    obj.dtc = sliding.total_nucleations - obj.total_nucleations
                    obj.side = side 
                    graph.add_node( move,
                                    obj = obj, 
                                    pairs = obj.total_nucleations,
                                    state = obj.state,
                                    dpxdist = obj.dpxdist,
                                    dtc = obj.dtc,
                                    side = obj.side)
                    graph.add_edge(structure, move) 
                    graph.add_edge(move, structure)
                    yield list(leafs(move, target, graph))
            #perform leafing for each backfray and connect everything in a unique graph
            subSG = nx.DiGraph()
            for bf in sliding.backfray:
                bfG = nx.DiGraph()
                list(leafs(bf.structure, sliding.structure, bfG))
                subSG.update(bfG)
            return subSG
                
        split = self.chamber.split_slidings()
        for l, r in zip(split[0], split[1]):
            lsubSG = backfraygraph(l, 'l')
            rsubSG = backfraygraph(r, 'r')
            for nl, nr in zip(lsubSG.nodes(), rsubSG.nodes()):
                lsubSG.nodes[nl]['G'] = lsubSG.nodes[nl]['obj'].structureG()
                rsubSG.nodes[nr]['G'] = rsubSG.nodes[nr]['obj'].structureG()
            for graph in [lsubSG, rsubSG]:
                most = list(self.filternodes('G', min, graph).nodes())[0]
                graph.nodes[most]['mostable'] = True
            self.SG.update(lsubSG)
            self.SG.update(rsubSG)
        
        for node in self.SG.nodes():
            self.SG.nodes[node]['obj'].structureG()

########CONNECT NUCLEATIONS TO SIMPLEX
        nucleations = self.filternodes('pairs', lambda x: x == self.min_nucleation, self.SG)
        for node in list(nucleations.nodes()):
            state = 'off_nucleation'
            self.SG.nodes[node]['state'] = state
            dgss = 0 #reference simplex state
            dgnuc = self.SG.nodes[node]['obj'].G #nucleation free energy
            # compute forward and backward rates
            fwd, bwd = self.nmethod(state, dgss, dgnuc)
            # normalize the forward nucleation rate
            # osl = self.ownsliding(self.SG.nodes[node], self.SG) 
            # bpt = self.SG.nodes[osl]['obj'].total_nucleations 
            nucnorm = (self.s1.length + self.s2.length - 1) * (self.chamber.duplex.total_nucleations - self.SG.nodes[node]['dpxdist'])
            fwd = fwd / nucnorm
            bwd = bwd / nucnorm
            self.SG.add_edge(self.simplex, node, k = fwd, state = state)
            self.SG.add_edge(node, self.simplex, k = bwd, state = state)
            next = list(self.SG.neighbors(node))
            next.remove(self.simplex)
            for e in next:
                dge = self.SG.nodes[e]['obj'].G
                # compute forward and backward rates
                fwd, bwd = self.zmethod('zipping', dgnuc, dge) 
                self.SG.add_edge(node, e, k = fwd, state = 'zipping')
                self.SG.add_edge(e, node, k = bwd, state = 'zipping')
        
########MAKE BACKFRAY CONNECTIONS
        mbare = list(set(self.SG.nodes()) - set(nucleations.nodes()))
        for e1, e2 in nx.subgraph_view(self.SG, lambda x: x in mbare).to_undirected().edges(): 
            state = 'zipping'
            dg1 = self.SG.nodes[e1]['obj'].G
            dg2 = self.SG.nodes[e2]['obj'].G
            fwd, bwd = self.zmethod(state, dg1, dg2)
            self.SG[e1][e2]['k'] = fwd
            self.SG[e2][e1]['k'] = bwd
            self.SG[e1][e2]['state'] = state
            self.SG[e2][e1]['state'] = state
        
########MAKE SLIDING CONNECTIONS
        left = self.filternodes('side', lambda x: x == 'l', self.SG)
        right = self.filternodes('side', lambda x: x == 'r', self.SG)
        lcompl = self.filternodes('mostable', lambda x: x == True, left)
        rcompl = self.filternodes('mostable', lambda x: x == True, right)
        for l, r in combinations([*list(lcompl.nodes.data()), *list(rcompl.nodes.data())], 2):
            state = 'sliding'
            self.SG.nodes[l[0]]['state'] = state
            self.SG.nodes[r[0]]['state'] = state
            #CONNECT SLIDINGS TO DUPLEX
            if l[1]['dpxdist'] <= self.sliding_cutoff:
                dgduplex = self.chamber.duplex.G
                dgsliding = self.SG.nodes[l[0]]['obj'].G
                fwd, bwd = self.smethod(state, dgsliding, dgduplex)
                fwd = fwd/l[1]['dpxdist']
                bwd = bwd/l[1]['dpxdist']
                self.SG.add_edge(l[0], self.duplex, k = fwd, state = state)
                self.SG.add_edge(self.duplex, l[0], k = bwd, state = state)
            if r[1]['dpxdist'] <= self.sliding_cutoff:
                dgduplex = self.chamber.duplex.G
                dgsliding = self.SG.nodes[r[0]]['obj'].G
                fwd, bwd = self.smethod(state, dgsliding, dgduplex)
                fwd = fwd/r[1]['dpxdist']
                bwd = bwd/r[1]['dpxdist']
                self.SG.add_edge(r[0], self.duplex, k = fwd, state = state)
                self.SG.add_edge(self.duplex, r[0], k = bwd, state = state)
            #CONNECT SLIDINGS BETWEEN THEMSELVES
            if l[1]['side'] == r[1]['side']:
                distance = abs(l[1]['dpxdist'] - r[1]['dpxdist'])
                if distance <= self.sliding_cutoff:
                    dgl = self.SG.nodes[l[0]]['obj'].G
                    dgr = self.SG.nodes[r[0]]['obj'].G
                    fwd, bwd = self.smethod(state, dgl, dgr)
                    fwd/distance
                    bwd/distance
                    self.SG.add_edge(l[0], r[0], k = fwd, state = state)
                    self.SG.add_edge(r[0], l[0], k = bwd, state = state)
            else:
                distance = abs(l[1]['dpxdist'] + r[1]['dpxdist'])
                if distance <= self.sliding_cutoff:
                    dgl = self.SG.nodes[l[0]]['obj'].G
                    dgr = self.SG.nodes[r[0]]['obj'].G
                    fwd, bwd = self.smethod(state, dgl, dgr)
                    fwd/distance
                    bwd/distance
                    self.SG.add_edge(l[0], r[0], k = fwd, state = state)
                    self.SG.add_edge(r[0], l[0], k = bwd, state = state)

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
    
    def ownsliding(self, node, graph):
        """ special filtering method to get the goal sliding associated with an off_nucleation"""
        dpxdist = node['dpxdist']
        side    = node['side']
        def filternode(node):
            try: return graph.nodes[node]['dpxdist'] == dpxdist and graph.nodes[node]['side'] == side and graph.nodes[node]['state'] == 'sliding'
            except KeyError: pass     
        sliding = list(nx.subgraph_view(graph, filter_node=filternode).nodes())[0]
        return sliding
        
    def slidingcondition(self, slide0, slide1, duplexation = False):
        # if duplexation == True: 
        #     if slide0.dpxdist <= self.sliding_cutoff: return True 
        #     else: return False 
        if slide0.dpxdist - slide1.dpxdist < self.sliding_cutoff: return True 
        else: return False  

    def get_traps(self):
        sab='+'.join([self.s1.sequence,self.s2.sequence])
        def addpar(string, i, char):
            return string[:i] + char + string[i+1:]
        wc = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        ss = self.chamber.singlestranded.structure.split('+')[0]
        for i, e1 in enumerate(self.s1.sequence):
            for j, e2 in enumerate(self.s2.sequence):
                if e1 == wc[e2]:
                    if i+j != len(ss)-1:
                        dpxdist = abs(i+j - len(ss) + 1)
                        l = addpar(ss, i, '(')
                        r = addpar(ss, j, ')')
                        trap = '+'.join([l,r])
                        obj = Complex(self.model, self.s1, self.s2, state='off_nucleation', structure = trap, dpxdist=dpxdist)
                        print(sab)
                        print(trap)
                        self.DG.add_node(   trap,
                                            obj = obj, 
                                            pairs = obj.total_nucleations,
                                            state = obj.state,
                                            dpxdist = obj.dpxdist,
                                            fre = obj.G)
                        dgss = 0
                        dgtrap = obj.G
                        fwd, bwd = self.nmethod('off_nucleation', dgss, dgtrap)
                        nucnorm = (self.s1.length + self.s2.length - 1) 
                        fwd = fwd/nucnorm
                        bwd = bwd/nucnorm
                        self.DG.add_edge(self.simplex, trap, k = fwd)
                        self.DG.add_edge(trap, self.simplex, k = bwd)

    def get_neighbor_zippings(self, structure, onlyup = False, onlydown = False):

        left, right = structure.split('+')

        def get_i(lst):
            indices = []
            for i, el in enumerate(lst):
                if i > 0 and el != lst[i-1]:
                    indices.append(i)
            return indices

        def update_structure(string, character: str):
            indices = get_i(string)
            indices_inv = get_i(string[::-1])
            updated = string
            for index in indices:
                updated = updated[:index-1] + character + updated[index:]
            updated_inv = updated[::-1]
            for index in indices_inv:
                updated_inv = updated_inv[:index-1] + character + updated_inv[index:]
            return updated_inv[::-1]

        up = '+'.join([update_structure(left,'('),update_structure(right,')')])
        down = '+'.join([update_structure(left,'.'),update_structure(right,'.')])

        if onlyup == True:
            return up
        elif onlydown == True:
            return down
        else: return up, down
    
    def save_graph(self, PATH):
        import os 
        #convert node object to string of object type
        for n in self.DG.nodes.data():
            n[1]['obj'] = str(type(n[1]['obj']))
        try: os.makedirs(PATH)
        except FileExistsError: pass 
        nx.write_gexf(self.DG,f'{PATH}/{self.s1.sequence}_graph_K.gexf')

    @property
    def nodes(self):
        return self.Graph.nodes()
    
    def nodata(self, *attributes):
        if attributes != None: 
            return list(self.Graph.nodes.data(*attributes))
        else: return list(self.Graph.nodes.data())


########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################


class Kinetics(object):
    def __init__(self, model: Model, chamber: Chamber):
        
        self.model = model
        self.T = model.kelvin
        self.space = model.space_dimensionality
        self.s1 = chamber.s1
        self.s2 = chamber.s2 
        
        self.phys = {'R(kcal/molK)': 1.987e-3,
                     'h_planck':     6.62607015e-34,
                     'k_boltz':      1.380649e-23,
                     'k_boltz_cm':   1.380649e-19,
                     'Na':           6.023e23,
                     'gamma':        .57722,}
        
        self.hydroparams = {'H20_viscosity': {'value':0.8701e-5, 'units':'kg/(cm*s)'}}
        self.viscosity = self.hydroparams['H20_viscosity']['value']

        #Lipid diffusion constants
        self.lipid_diffusion = {'Units':'cm^2/s',
                                'Filippov04': 9.32e-8,
                                'Arnott08':   4.5e-8}

        #TODO DOUBLE CHECK THESE VALUES 
        self.duplex_geometry = {'units':'cm',
                                'basepairdistance': 0.34e-7, 
                                'persistence_length': 100*0.34e-7,
                                'radius': 2e-7} #CHECK THE RADIUS 

        #TODO DOUBLE CHECK THESE VALUES 
        self.simplex_geometry = {   'units':'cm',
                                    'basepairdistance': 0.676e-7,
                                    'persistence_length': 2.223e-7,
                                    'cylinder_radius': 1e-7}

        ### Default compute relevant rates
        self.diffusionlimited()
        self.geometric_rate()

    
        """ ------------- SLIDING FORWARD RATES -----------------"""

    def set_slidingrate(self, sliding): #NOTE NEED TO EXPAND THIS TO TRY ALL RANGES
        self.slidingrate = sliding
    

    """ ------------- ZIPPING FORWARD RATES -----------------"""

    def set_zippingrate(self, setrate):
        """ first approximation is to take zipping
            equal to diffusion limited collision rate """
        if setrate == 'diffusionlimited': self.zippingrate = self.dlrate
        else: self.zippingrate = setrate
    
    
    def frates(self, kind):
        correspondence = {'f_on_core':      self.georate,
                          'f_zipping':      self.zippingrate,
                          'f_zipping-end':  self.zippingrate,
                          'f_off_core':     self.georate,
                          'f_backfray':     self.zippingrate,
                          'f_backfray-end': self.zippingrate,
                          'f_sliding':      self.slidingrate,
                          'f_sliding-end':  self.slidingrate,
                        }
        return correspondence[kind]
    
    def brates(self, kind, dg):
        fkind = 'f'+kind[1:]
        kf = self.frates(fkind)
        correspondence = {'b_on_core':      self.k_back(kf, dg),
                          'b_zipping':      self.k_back(kf, dg),
                          'b_zipping-end':  self.k_back(kf, dg),
                          'b_off_core':     self.k_back(kf, dg),
                          'b_backfray':     self.k_back(kf, dg),
                          'b_backfray-end': self.k_back(kf, dg),
                          'b_sliding':      self.k_back(kf, dg),
                          'b_sliding-end':  self.k_back(kf, dg)
                        }
        return correspondence[kind]
    
    def unif_scaling(self, nucleations):
        return 1/nucleations

    def z_scaling(self, nucleations):
        """
        Returns:
        - the nucleation partition function as sum over nodes of e^dg_node/KbT
        - list of boltzmann weights per node (use a dictionary)
        Then boltzmann weighting will just be indexing the given node
        from the dictionary and dividing the value by Z 
        """
        pass
    

    def kawasaki(self, kind, dgi, dgj):
        ratesdict = {'zipping':     self.zippingrate,
                     'backfray':    self.zippingrate,
                     'duplex':      self.zippingrate,
                     'sliding':     self.slidingrate,
                     'on_nucleation':   self.georate,
                     'off_nucleation':  self.georate}
        ka = ratesdict[kind]
        #kf
        kij = ka*np.exp(-(dgj-dgi)/(2*self.phys['R(kcal/molK)']*(self.T)))
        #kb
        kji = ka*np.exp(-(dgi-dgj)/(2*self.phys['R(kcal/molK)']*(self.T)))
        return kij, kji
    
    def metropolis(self, rate, dgi, dgj):
        ratesdict = {'zipping':     self.zippingrate,
                     'backfray':    self.zippingrate,
                     'duplex':      self.zippingrate,
                     'sliding':     self.slidingrate,
                     'on_nucleation':   self.georate,
                     'off_nucleation':  self.georate}
        ka = ratesdict[rate]
        if dgi > dgj:
            kij = ka
            kji = ka * np.exp(-(dgi-dgj)/(self.phys['R(kcal/molK)']*(self.T)))
        elif dgi < dgj:
            kij = ka * np.exp(-(dgj-dgi)/(self.phys['R(kcal/molK)']*(self.T)))
            kji = ka 
        else: 
            kij = ka
            kji = ka 
        return kij, kji
            

    ##################################################
    #################### General #####################
    ##################################################

    def einsmol_spherical(self, radius):
        """ returns einstein smoluchowski diffusivity for 
            spherical particles in (cm^2)/s units """
        return self.phys['k_boltz_cm']*self.T/(6*np.pi*self.viscosity*radius)

    def diffusionlimited(self, kind='ppi'):
        """ Smoluchowski 1916 classical formula """
        if kind == 'ppi':  
            self.ppi_diffusivities()
            dc = 1e-3 #dimensional correction from cubic cm to cubic dm to get 1/(M*s) kinetic rates
            self.dlrate = self.phys['Na']*4*np.pi*(self.pD1+self.pD2)*(self.gr1+self.gr2) * dc
            return self.dlrate
        elif kind == 'vanilla':
            self.vanilla_diffusivities()
            dc = 1e-3 #dimensional correction from cubic cm to cubic dm to get 1/(M*s) kinetic rates
            self.dlrate = self.phys['Na']*4*np.pi*(self.vD1+self.vD2)*(self.size1+self.size2) * dc
            return self.dlrate
        else: raise ValueError(f'{kind} not implemented')

    def geometric_rate(self):
        self.georate = (np.power((self.model.theta/360),2))*(np.power((self.model.phi/360),2)) * self.dlrate
        return self.georate 

    def closedconfscaling(self, p_circular):
        """ p_circular => probability of circularization 
            circularization will impede nucleation because ...
            (just make a picture of it in your mind for now) """
        self.georate = self.georate * (1 - p_circular)
        #TODO implement methods to calculate p_circular from 
        #end-to-end probability distributions from polymer physics 



    """ ------------- NUCLEATION FORWARD RATES -----------------"""

    ##################################################
    #################### Vanilla #####################
    ##################################################
   
    def ss_strands_size(self):
        self.size1 = self.s1.length*self.simplex_geometry['basepairdistance']
        self.size2 = self.s2.length*self.simplex_geometry['basepairdistance']
        # TODO DIMENSIONAL ANALYSIS HERE (RIGHT NOW IT IS WRONG)
   
    def vanilla_diffusivities(self):
        self.ss_strands_size()
        self.vD1 = self.einsmol_spherical(self.size1)
        self.vD2 = self.einsmol_spherical(self.size2)


    #########################################
    ####### Polymer Physics Informed ########
    #########################################
    
    def gyradiuses(self):
        lambda_0 = (self.simplex_geometry['persistence_length']*self.simplex_geometry['basepairdistance'])/3
        self.gr1 = np.sqrt(self.s1.length*lambda_0)
        self.gr2 = np.sqrt(self.s2.length*lambda_0)
    
    def ppi_diffusivities(self):
        self.gyradiuses()
        self.pD1 = self.einsmol_spherical(self.gr1)
        self.pD2 = self.einsmol_spherical(self.gr2)
    

  
    #TODO###################################################
    ######## Chain wiggling pdfs to implement later ########
    def px_realchain(x):
        return 0.278*(x**0.28)*np.exp(-1.206*(x**2.43))
    def px_idealchain(x):
        return ((3/2*np.pi)**(3/2))*np.exp(-1.5*(x**2))
    
    """ Check dnapolymer.py inside ./oldsequentialcode/polymer 
        for all the other polymer physics related functions """
        

    """ ------------- GENERAL BACKWARD RATES -----------------"""

    def k_equilibrium(self, free_energy):
        return np.exp(-(free_energy/(self.phys['R(kcal/molK)']*(self.T))))   # (Kcal/mol)/((Kcal/mol*K)*K) -> adimensional
    
    def k_back(self, forward, free_energy, geo = 'cylinder', p_circular = None):
        """ LOOK OUT: default angle steric values don't influence 
            the backward rates (since they are yelding a factor of 1)"""
        ke = self.k_equilibrium(free_energy)
        self.diffusionlimited()
        self.geometric_rate()
        if geo == 'cylinder': return self.georate / ke 
        if geo == 'chain:': 
            if 0 <= p_circular <= 1 : return self.closedconfscaling(p_circular) / ke
            else: raise ValueError("need to input a 'p_circular' value in the [0,1] interval")
        else: return forward / ke


    #############################################
    ########## General helper methods ###########
    #############################################
    


class FatalError(Exception):
    def __init__(self, message):
        super().__init__(message)

