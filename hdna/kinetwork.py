import networkx as nx
import numpy as np 
import pandas as pd 

from itertools import pairwise, combinations, tee

from .complex import Complex
from .model import Model
from .strand import Strand
 
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

        self.nmethod = self.kinetics.kawasaki
        self.zmethod = self.kinetics.kawasaki
        self.smethod = self.kinetics.metropolis

        self.normalizeback = False
        self.nucnorm = np.power((self.s1.length + self.s2.length - self.model.min_nucleation + 1),2)

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
        ss = self.simplex 
        num = min(self.s1.length, self.s2.length)
        for n in range(num):
            for i, e1 in enumerate(self.nwise(self.s1.sequence, n)):
                for j, e2 in enumerate(self.nwise(self.s2.sequence, n)):
                    e1 = self.u(*e1)
                    e2 = self.u(*e2)
                    if e1 == self.wc(e2[::-1]):
                         # --> Condition for checking that the nucleation is off register 
                        if verbose: print(self.sab(self.s1.sequence, self.s2.sequence))
                        spacei = ' '*(i)
                        spacej = ' '*(j - i + len(ss) - n + 1)
                        if verbose: print(spacei+spacej.join([e1,e2]))
                        dpxdist = len(ss) - n - i - j
                        l = self.addpar(ss, i, n, '(')
                        r = self.addpar(ss, j, n, ')')
                        trap = self.sab(l,r)
                        if i+j != len(ss)-n:
                            state = 'off_nucleation' if n == 1 else 'backfray' 
                        else: state = 'on_nucleation' if n == 1 else 'zipping'
                        obj = Complex(self.model, self.s1, self.s2, state=state, structure = trap, dpxdist=dpxdist)
                        self.DG.add_node(   trap,
                                            obj = obj, 
                                            pairs = obj.total_nucleations,
                                            state = obj.state,
                                            dpxdist = obj.dpxdist,
                                            tdx =(i,j),
                                            fre = obj.structureG())
                        dgss = 0
                        dgtrap = obj.G
                        fwd, bwd = self.nmethod('off_nucleation', dgss, dgtrap)
                        fwd = fwd / self.nucnorm
                        if self.normalizeback: bwd = bwd / self.nucnorm
                        if n == 1:
                            self.DG.add_edge(self.simplex, trap, k = fwd, state = 'off_nucleation')
                            self.DG.add_edge(trap, self.simplex, k = bwd, state = 'off_nucleation')
                        elif n > 1:
                            self.sldbranches.append(dpxdist)
                            f1 = self.filternodes('dpxdist', lambda x: x == obj.dpxdist, self.DG)
                            f2 = self.filternodes('pairs', lambda x: x == n-1, f1)
                            f3 = self.filternodes('tdx', lambda x: (i-1 <= x[0] <= i+1) and (j-1 <= x[1] <= j+1), f2)
                            for node in f3.nodes():
                                if verbose: print(node, trap)
                                dgnode = self.DG.nodes[node]['fre']
                                fwd, bwd = self.zmethod('zipping', dgnode, dgtrap)
                                self.DG.add_edge(node, trap, k = fwd, state = 'backfray')
                                self.DG.add_edge(trap, node, k = bwd, state = 'backfray') 
                            if n == num-1: 
                                dgduplex = self.DG.nodes[self.duplex]['fre']
                                fwd, bwd = self.zmethod('zipping', dgtrap, dgduplex)
                                self.DG.add_edge(trap, self.duplex, k = fwd, state = 'backfray')
                                self.DG.add_edge(self.duplex, trap, k = bwd, state = 'backfray') 
                        
        #get unique ids for all branches except branch 0 corresponding to on register nucleation
        self.sldbranches = set(self.sldbranches)
        self.sldbranches.remove(0)

##########################################################################
##########################################################################

    def connect_slidings(self):
        for branch in self.sldbranches:
            leaf = self.filternodes('dpxdist', lambda x: x == branch, self.DG)
            mostable = list(self.filternodes('fre', min, leaf).nodes())[0]
            self.DG.nodes[mostable]['state'] = 'sliding'
            dgsliding = self.DG.nodes[mostable]['fre']
            dgduplex = self.DG.nodes[self.duplex]['fre']
            fwd, bwd = self.smethod('sliding', dgsliding, dgduplex)
            fwd = fwd / abs(branch)
            bwd = bwd / abs(branch)
            self.DG.add_edge(mostable, self.duplex, k = fwd, state = 'sliding')
            self.DG.add_edge(self.duplex, mostable, k = bwd, state = 'sliding')

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
    
##########################################################################

    def save_graph(self, PATH):
        import os 
        #convert node object to string of object type
        for n in self.DG.nodes.data():
            n[1]['obj'] = str(type(n[1]['obj']))
            try: n[1]['tdx'] = str(type(n[1]['tdx']))
            except KeyError: pass
        try: os.makedirs(PATH)
        except FileExistsError: pass 
        nx.write_gexf(self.DG,f'{PATH}/{self.s1.sequence}_graph_K.gexf')

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
        return self.Graph.nodes()
    
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


########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################


class Kinetics(object):
    def __init__(self, model: Model, s1: Strand, s2: Strand):
        
        self.model = model
        self.T = model.kelvin
        self.space = model.space_dimensionality
        self.s1 = s1
        self.s2 = s2 
        
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

