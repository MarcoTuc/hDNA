import networkx as nx
import numpy as np 
import pandas as pd 

from itertools import pairwise, chain, tee

from .chamber import Chamber
from .model import Model, Geometry
from .strand import Strand
 
class Kinetwork(object):

    """Given a chamber it will generate the corresponding
    Kinetic network to be later passed to the simulator"""

    def __init__(self, model: Model, s1: Strand, s2: Strand, geometry: Geometry, clean=False):

        self.model = model 
        self.s1 = s1        #53
        self.s2 = s2        #35
        self.min_nucleation = self.model.min_nucleation
        self.sliding_cutoff = self.model.sliding_cutoff
        self.chamber = Chamber(self.model, self.s1, self.s2)
        self.kinetics = Kinetics(self.model, self.chamber, geometry)
        self.kinetics.set_slidingrate(self.model.sliding)
        self.kinetics.set_zippingrate(self.model.zipping)
        self.Graph = nx.DiGraph()
        self.add_nodes()
        self.relabel_nodes()

        if not clean:           
            self.add_reactions()
            self.clean_duplicates_relabel()

            #get an overview of the network 
        self.possible_states = [ None,
                                'singlestranded',
                                'duplex',
                                'zipping',
                                'on_nucleation',
                                'off_nucleation', 
                                'backfray',
                                'sliding' ]   

        countslist = [list(dict(self.Graph.nodes.data('state')).values()).count(state) for state in self.possible_states]
        self.overview = dict(zip(self.possible_states,countslist))


    def add_nodes(self):

        # Make a list of notes object from chamber
        nodes = [self.chamber.singlestranded]
        for s in self.chamber.slidings:
            nodes.append(s)
            for bf in s.backfray:
                nodes.append(bf)
                for z in bf.zipping:
                    nodes.append(z)
        for s in self.chamber.oncores:
            nodes.append(s)
            for z in s.zipping:
                nodes.append(z)
        nodes.append(self.chamber.duplex)
        
        # Chamber makes (lots of) duplicate objects when 
        # parsing structures, here I remove that shit 
        DF = pd.DataFrame(nodes)
        DF.columns = ['object']
        DF['structure'] = [e.structure for e in DF['object']]
        DF = DF.drop_duplicates(subset='structure', keep='first')
        DF = DF.reset_index()
        DF = DF.drop('index', axis=1)
        
        for i, row in DF.iterrows():
            try: statecheck = row['object'].state
            except: print(row['object'], row['object'].structure)
            self.Graph.add_node(
                row['object'],
                obj = row['object'],
                structure = row['structure'],
                state = statecheck,
                pairs = row['object'].total_nucleations,
                dpxdist = row['object'].dpxdist
            )

    def relabel_nodes(self):
        struct = [e['structure'] for e in self.Graph._node.values()]
        labels = self.Graph._node.keys()
        diz = dict(zip(labels, struct))
        Relabeled = nx.relabel_nodes(self.Graph, diz)
        self.Graph = Relabeled


    def add_reactions(self, verbose=False):

        def pairwisend(iterable):
            from itertools import zip_longest
            a, b = tee(iterable)
            next(b, None)
            return zip_longest(a, b)
        
        simplex = self.chamber.singlestranded.structure
        duplex  = self.chamber.duplex.structure
        
        df = pd.DataFrame(self.nodata('state'), columns=['structure','state'])
        
        native = df[df['state'] == 'on_nucleation']
        for n in native['structure']:
            cplx = self.nodes[n]['obj']
            self.Graph.add_edge(simplex, 
                                n, 
                                kind = 'f_on_core', 
                                k = self.kinetics.frates('f_on_core'))
            self.Graph.add_edge(n, 
                                simplex, 
                                kind = 'b_on_core', 
                                k = self.kinetics.brates('b_on_core', 
                                                        cplx.G))
            self.Graph.add_edge(n, 
                                cplx.zipping[0].structure, 
                                kind = 'f_zipping',
                                k = self.kinetics.frates('f_zipping'))
            self.Graph.add_edge(cplx.zipping[0].structure, 
                                n, 
                                kind = 'b_zipping',
                                k = self.kinetics.brates('b_zipping',
                                                        cplx.zipping[0].G - cplx.G))
            for z1, z2 in pairwise(cplx.zipping):
                self.Graph.add_edge(z1.structure, 
                                    z2.structure, 
                                    kind = 'f_zipping',
                                    k = self.kinetics.frates('f_zipping'))
                self.Graph.add_edge(z2.structure, 
                                    z1.structure, 
                                    kind = 'b_zipping',
                                    k = self.kinetics.brates('b_zipping',
                                                            z2.G - z1.G))
            self.Graph.add_edge(z2.structure, 
                                duplex, 
                                kind = 'f_zipping-end',
                                k = self.kinetics.frates('f_zipping-end'))
            self.Graph.add_edge(duplex, 
                                z2.structure, 
                                kind = 'b_zipping-end',
                                k = self.kinetics.brates('b_zipping-end',
                                                        self.chamber.duplex.G - z2.G))
        
        nonnative = df[df['state'] == 'sliding']      
         
        for n, nxt in pairwisend(nonnative['structure']):

            cplx = self.nodes[n]['obj']
            if nxt != None:
                cnxt = self.nodes[nxt]['obj']

            if not cplx.bfempty:
                #FORWARD
                self.Graph.add_edge(simplex, 
                                    cplx.backfray[0].structure, 
                                    kind = 'f_off_core',
                                    k = self.kinetics.frates('f_off_core'))
                if verbose: print(self.Graph.has_edge(simplex, cplx.backfray[0].structure))
                #BACKWARD
                self.Graph.add_edge(cplx.backfray[0].structure, 
                                    simplex, 
                                    kind = 'b_off_core',
                                    k = self.kinetics.brates('b_off_core',
                                                            cplx.backfray[0].G))
                if verbose: print(self.Graph.has_edge(cplx.backfray[0].structure, simplex))
                for bf in cplx.backfray:
                    if verbose: 
                        print('NEWBF:')
                        print('o', n)   
                        print('b', bf.structure)
                    #FORWARD
                    self.Graph.add_edge(simplex, 
                                        bf.structure, 
                                        kind = 'f_off_core',
                                        k = self.kinetics.frates('f_off_core'))
                    if verbose: print(self.Graph.has_edge(simplex, bf.structure))
                    #BACKWARD
                    self.Graph.add_edge(bf.structure, 
                                        simplex, 
                                        kind = 'b_off_core',
                                        k = self.kinetics.brates('b_off_core', 
                                                                bf.G))
                    if verbose: 
                        print(self.Graph.has_edge(bf.structure,simplex))
                        print(len(bf.zipping))
                    if len(bf.zipping) > 1:
                        if verbose: 
                            print('n',n)
                            print('b',bf.structure)
                        for z in bf.zipping:
                            if verbose: print('z',z.structure)
                        #FORWARD
                        self.Graph.add_edge(bf.structure, 
                                            bf.zipping[0].structure, 
                                            kind = 'f_backfray',
                                            k = self.kinetics.frates('f_backfray'))
                        #BACKWARD
                        self.Graph.add_edge(bf.zipping[0].structure, 
                                            bf.structure, 
                                            kind = 'b_backfray',
                                            k = self.kinetics.brates('b_backfray',
                                                                    bf.zipping[0].G - bf.G))
                        if verbose: print('im here in >1')
                        for z1, z2 in pairwise(bf.zipping): 
                            if verbose: 
                                print('z', z1.structure)
                                print('z', z2.structure)
                            #FORWARD
                            self.Graph.add_edge(z1.structure, 
                                                z2.structure, 
                                                kind = 'f_backfray',
                                                k = self.kinetics.frates('f_backfray'))
                            #BACKWARD
                            self.Graph.add_edge(z2.structure, 
                                                z1.structure, 
                                                kind = 'b_backfray',
                                                k = self.kinetics.brates('b_backfray', 
                                                                        z2.G - z1.G))
                        #FORWARD
                        self.Graph.add_edge(bf.zipping[-1].structure, 
                                            n, 
                                            kind = 'f_backfray-end',
                                            k = self.kinetics.frates('f_backfray-end'))
                        #BACKWARD
                        self.Graph.add_edge(n, 
                                            bf.zipping[-1].structure, 
                                            kind = 'b_backfray-end',
                                            k = self.kinetics.brates('b_backfray-end',
                                                                    cplx.G - bf.zipping[-1].G))
                    elif len(bf.zipping) == 1:
                        if verbose: print('im here in 1')
                        #BF0TOITSZIPP
                        #FORWARD
                        self.Graph.add_edge(bf.structure, 
                                            bf.zipping[0].structure, 
                                            kind = 'f_backfray-end',
                                            k = self.kinetics.frates('f_backfray-end'))
                        #BACKWARD
                        self.Graph.add_edge(bf.zipping[0].structure, 
                                            bf.structure,
                                            kind = 'b_backfray-end',
                                            k = self.kinetics.brates('b_backfray-end',
                                                                    cplx.G - bf.zipping[0].G))    
                        #ZIPPTOSTRUCT
                        #FORWARD
                        self.Graph.add_edge(bf.zipping[0].structure, 
                                            n, 
                                            kind = 'f_backfray-end',
                                            k = self.kinetics.frates('f_backfray-end'))
                        #BACKWARD
                        self.Graph.add_edge(n, 
                                            bf.zipping[0].structure, 
                                            kind = 'b_backfray-end',
                                            k = self.kinetics.brates('b_backfray-end',
                                                                    cplx.G - bf.zipping[0].G))    
                    elif len(bf.zipping) == 0:
                        if n == bf.structure:
                            if verbose: print('im here in 0')
                        #FORWARD
                        self.Graph.add_edge(bf.structure, 
                                            n, 
                                            kind = 'f_backfray-end',
                                            k = self.kinetics.frates('f_backfray-end'))
                        #BACKWARD
                        self.Graph.add_edge(n, 
                                            bf.structure, 
                                            kind = 'b_backfray-end',
                                            k = self.kinetics.brates('b_backfray-end',
                                                                    cplx.G - bf.G))
            else: 
                #FORWARD
                self.Graph.add_edge(simplex, 
                                    n, 
                                    kind = 'f_backfray-end',
                                    k = self.kinetics.frates('f_off_core'))
                #BACKWARD
                self.Graph.add_edge(n, 
                                    simplex, 
                                    kind = 'b_backfray-end',
                                    k = self.kinetics.brates('b_off_core',
                                                            cplx.G))
            if self.slidingcondition(cplx, self.chamber.duplex):
                #FORWARD
                self.Graph.add_edge(n, 
                                    duplex, 
                                    kind = 'f_sliding-end',
                                    k = self.kinetics.frates('f_sliding-end'))
                #BACKWARD
                self.Graph.add_edge(duplex, 
                                    n, 
                                    kind = 'b_sliding-end',
                                    k = self.kinetics.brates('b_sliding-end',
                                                            self.chamber.duplex.G - cplx.G))
                
            if self.slidingcondition(cplx, cnxt) and nxt != None:
                if cplx.dpxdist > cnxt.dpxdist:
                    #FORWARD
                    self.Graph.add_edge(n, 
                                        nxt, 
                                        kind = 'f_sliding',
                                        k = self.kinetics.frates('f_sliding'))
                    #BACKWARD
                    self.Graph.add_edge(nxt, 
                                        n, 
                                        kind = 'b_sliding',
                                        k = self.kinetics.brates('b_sliding',
                                                                cnxt.G - cplx.G))
                else:                 
                    self.Graph.add_edge(nxt, 
                                        n, 
                                        kind = 'f_sliding',
                                        k = self.kinetics.frates('f_sliding'))
                    self.Graph.add_edge(n, 
                                        nxt, 
                                        kind = 'b_sliding',
                                        k = self.kinetics.brates('b_sliding',
                                                                cplx.G - cnxt.G))
                
            if verbose: print('\n')



        
    def slidingcondition(self, slide0, slide1, duplexation = False):
        # if duplexation == True: 
        #     if slide0.dpxdist <= self.sliding_cutoff: return True 
        #     else: return False 
        if slide0.dpxdist - slide1.dpxdist < self.sliding_cutoff: return True 
        else: return False  


    def clean_duplicates_relabel(self):
        """ Unfortunately my code makes duplicates.
            Also unfortunately it was easier to just get rid of them
            after they are made than fix the entire code.
            Also the labels will be changed from object to structure.
            Objects are still accessible from nodes object property. """
        df = pd.DataFrame(self.Graph.nodes(data=True))
        labels = df[0]
        struct = [list(df[1])[i]['structure'] for i in range(len(list(df[1])))]
        diz = dict(zip(labels, struct))
        Relabeled = nx.relabel_nodes(self.Graph, diz)
        self.Graph = Relabeled

    def node_filter(self, property, attribute):
        return self.Graph.subgraph( 
        [n for n, attrdict in self.Graph.nodes.items() if attrdict [str(property)] == str(attribute)])

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
        for n in self.Graph.nodes.data():
            n[1]['obj'] = str(type(n[1]['obj']))
        try: os.makedirs(PATH)
        except FileExistsError: pass 
        nx.write_gexf(self.Graph,f'{PATH}/{self.s1.sequence}_graph_K.gexf')

    @property
    def nodes(self):
        return self.Graph.nodes()
    
    def nodata(self, *attributes):
        if attributes != None: 
            return list(self.Graph.nodes.data(*attributes))
        else: return list(self.Graph.nodes.data())

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################



class Kinetics(object):
    def __init__(self, model: Model, chamber: Chamber, geometry: Geometry):
        
        self.model = model
        self.geometry = geometry
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
        self.set_slidingrate()
        self.set_zippingrate()
    
    
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
        self.georate = (np.power((self.geometry.theta/360),2))*(np.power((self.geometry.phi/360),2)) * self.dlrate
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
        

    """ ------------- SLIDING FORWARD RATES -----------------"""

    def set_slidingrate(self, sliding = 2e7): #NOTE NEED TO EXPAND THIS TO TRY ALL RANGES
        self.slidingrate = sliding
    

    """ ------------- ZIPPING FORWARD RATES -----------------"""

    def set_zippingrate(self, setrate=None):
        """ first approximation is to take zipping
            equal to diffusion limited collision rate """
        if setrate == None: self.zippingrate = self.dlrate
        elif type(setrate) == float: self.zippingrate = setrate


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
    
    # def rates(self, kind):
    #     correspondence = {'sliding':        self.slidingrate,
    #                       'sliding-end':    self.slidingrate,
    #                       'off_zipping':    self.zippingrate,
    #                       'zipping':        self.zippingrate,
    #                       'zipping-end':    self.zippingrate}
    #     return correspondence[kind]


class FatalError(Exception):
    def __init__(self, message):
        super().__init__(message)


#################################################################################################################
###########################################À GRAVEYARD ##########################################################
################################ MAY THESE FUNCTIONS REST IN PEACE ##############################################
################################################################################################################# 
  


                # self.Graph.add_node(
        #     self.chamber.singlestranded, 
        #     structure = self.chamber.singlestranded.structure, 
        #     state = 'singlestranded', 
        #     pairs = 0)
        # for s in self.chamber.slidings:
        #     check = len(s.backfray) > 1
        #     self.Graph.add_node(
        #         s, 
        #         structure = s.structure, 
        #         state = 'off_register', 
        #         pairs = int(s.total_nucleations), 
        #         dpxdist = s.dpxdist)
        #     if check:
        #         for bf in s.backfray: #first backfray is the offcored sliding 
        #             self.Graph.add_node(
        #                 bf,
        #                 structure = bf.structure,
        #                 state = 'backfray',
        #                 pairs = int(bf.total_nucleations))
        # for s in self.chamber.oncores:
        #     self.Graph.add_node(
        #         s,
        #         structure = s.structure,
        #         state = 'on_register',
        #         pairs = int(s.total_nucleations))
        # for s in self.chamber.oncores:
        #     for z in s.zipping[1:]: #first zipping value is the oncore itself so flush it away by indexing
        #         self.Graph.add_node(
        #             z, 
        #             structure = z.structure, 
        #             state = 'zipping', 
        #             pairs = int(z.total_nucleations))

        # self.Graph.add_node(
        #     self.chamber.duplex, 
        #     structure = self.chamber.duplex.structure, 
        #     state = 'duplex', 
        #     pairs = self.chamber.duplex.total_nucleations)
        # """   note that duplicated nodes will not be added to the resulting graph, hence I can run through 
        #       all zipping states and be sure that I will only get one note per each zipping that is common
        #       from different starting native nucleations. """


    #ADD RATES KINETWORK
        # ON  = self.node_filter('state','on_register')
        # D = self.node_filter('state', 'duplex')
        # duplex = list(D.nodes())[0]

        # for on, data in ON.nodes.items():
        #     self.Graph.add_edge(self.chamber.singlestranded, on, kind = 'on_nucleation') #ADDPROPERTY
        #     self.Graph.add_edge(on, on.zipping[1], kind = 'zipping')
        #     for i, z in enumerate(on.zipping[2:], start=2):
        #         self.Graph.add_edge(z, on.zipping[i-1], kind = 'zipping')
        #     self.Graph.add_edge(on.zipping[-1], self.chamber.duplex, kind = 'zipping-end')

        # L, R = self.chamber.split_slidings()
        # for l1, l2 in pairwise(L):
        #     print(l1.bfempty)
        #     print(l2.bfempty)
        #     if l1.bfempty:
        #         self.Graph.add_edge(self.chamber.singlestranded, l1, kind = 'off_nucleation')
        #     else:
        #         self.Graph.add_edge(self.chamber.singlestranded, l1.backfray[0], kind = 'off_nucleation')
        #         for bf1, bf2 in pairwise(l1.backfray[1:]):
        #             self.Graph.add_edge(bf1, bf2, kind = 'off_zipping') 
        #         if self.slidingcondition(l1.backfray[-1], duplex):
        #             self.Graph.add_edge(l1.backfray[-1], duplex, kind = 'off_zipping_end')

        #     if l2.bfempty: 
        #         self.Graph.add_edge(self.chamber.singlestranded, l2, kind = 'off_nucleation')
        #     else:
        #         self.Graph.add_edge(self.chamber.singlestranded, l2.backfray[0], kind = 'off_nucleation')
        #         for bf1, bf2 in pairwise(l2.backfray[1:]):
        #             self.Graph.add_edge(bf1, bf2, kind = 'off_zipping') 
        #         if self.slidingcondition(l2.backfray[-1], duplex):
        #             self.Graph.add_edge(l2.backfray[-1], duplex, kind = 'off_zipping_end')
            
        #     if self.slidingcondition(l1, l2):
        #         self.Graph.add_edge(l1, l2)
            


    # if False:
    #     L, R = self.chamber.split_slidings()
    #     if verbose:
    #         for l, r in zip(L, R):
    #             print(l.s1.sequence+'+'+l.s2.sequence)
    #             print(l.structure,'L')
    #             print('basepairs:',l.total_nucleations)
    #             print('dupdist:  ',l.dpxdist)
    #             print('\n')
    #             print(r.s1.sequence+'+'+r.s2.sequence)
    #             print(r.structure,'R')
    #             print('basepairs:',r.total_nucleations)
    #             print('dupdist:  ',r.dpxdist)
    #             print('\n') 
    #     for i, (left, right) in enumerate(zip(L, R)):
    #         lcheck = len(left.backfray) > 1
    #         rcheck = len(right.backfray) > 1
    #         if verbose:
    #             print('left: ',left.dpxdist)
    #             print('right:',right.dpxdist)
    #         if lcheck and rcheck: 
    #             self.Graph.add_edge(self.chamber.singlestranded, left.backfray[0], kind = 'off_nucleation')
    #             self.Graph.add_edge(self.chamber.singlestranded, right.backfray[0], kind = 'off_nucleation')
    #             for i, (bfl, bfr) in enumerate(zip(left.backfray[:-1], right.backfray[:-1])):
    #                 self.Graph.add_edge(bfl, left.backfray[i+1], kind = 'off_zipping')
    #                 self.Graph.add_edge(bfr, right.backfray[i+1], kind = 'off_zipping')
    #             if i > 0 and self.slidingcondition(left, L[i-1]) and self.slidingcondition(right, R[i-1]): 
    #                 self.Graph.add_edge(left, L[i-1], kind = 'sliding')
    #                 self.Graph.add_edge(right, R[i-1], kind = 'sliding')
    #         elif not lcheck and not rcheck:
    #             self.Graph.add_edge(self.chamber.singlestranded, left, kind = 'off_nucleation')
    #             self.Graph.add_edge(self.chamber.singlestranded, right, kind = 'off_nucleation')
    #             if i > 0 and self.slidingcondition(left, L[i-1]) and self.slidingcondition(right, R[i-1]): 
    #                 self.Graph.add_edge(left, L[i-1], kind = 'sliding')
    #                 self.Graph.add_edge(right, R[i-1], kind = 'sliding')
    #         elif not lcheck or not lcheck and not (lcheck and rcheck):
    #             raise FatalError('you messed up the slidings')
    #     try: 
    #         if self.slidingcondition(L[-1], None, duplexation=True): 
    #             self.Graph.add_edge(L[-1], list(D.nodes())[0], kind = 'sliding-end')
    #             if verbose: print('made sliding end connection')
    #     except IndexError: 
    #         if verbose: print('no left slidings as you can see from the empty list:', L)
    #         else: pass
    #     try: 
    #         if self.slidingcondition(R[-1], None, duplexation=True): 
    #             self.Graph.add_edge(R[-1], list(D.nodes())[0], kind = 'sliding-end')
    #             if verbose: print('made sliding end connection')
    #     except IndexError: 
    #         if verbose: print('no right slidings as you can see from the empty list:', R)
    #         else: pass