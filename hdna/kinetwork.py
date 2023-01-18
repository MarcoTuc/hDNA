import networkx as nx
import numpy as np 
import pandas as pd 

from hdna.chamber import Chamber
from hdna.model import Model, Geometry
from hdna.strand import Strand
 
class Kinetwork(object):

    """Given a chamber it will generate the corresponding
    Kinetic network to be later passed to the simulator"""

    def __init__(self, model: Model, s1: Strand, s2: Strand, minimum_nucleation):

        self.model = model 
        self.s1 = s1
        self.s2 = s2
        self.mincore = minimum_nucleation
        self.chamber = Chamber(self.model, self.s1, self.s2, self.mincore)
        

        self.Graph = nx.Graph()
        self.add_nodes()
        self.add_reactions()
        self.clean_duplicates()


    def add_nodes(self):

        self.Graph.add_node(self.chamber.singlestranded, object = self.chamber.singlestranded, structure = self.chamber.singlestranded.structure, state = 'singlestranded', pairs = 0)

        for s in self.chamber.offcores:
            self.Graph.add_node(s, object = s, structure = s.structure, state = 'off_register', pairs = int(s.total_nucleations))

        for s in self.chamber.oncores:
            self.Graph.add_node(s, object = s, structure = s.structure, state = 'on_register', pairs = int(s.total_nucleations))

        for s in self.chamber.oncores:
            for z in s.zipping[1:]: #first zipping value is the oncore itself so flush it away by indexing
                self.Graph.add_node(z, object = z, structure = z.structure, state = 'zipping', pairs = int(z.total_nucleations))

        self.Graph.add_node(self.chamber.duplex, object = self.chamber.duplex, structure = self.chamber.duplex.structure, state = 'duplex', pairs = self.chamber.duplex.total_nucleations)
        """   note that duplicated nodes will not be added to the resulting graph, hence I can run through 
              all zipping states and be sure that I will only get one note per each zipping that is common
              from different starting native nucleations. """


    def add_reactions(self, verbose=False):

        ON  = self.node_filter('state','on_register')
        D = self.node_filter('state', 'duplex')

        for on, data in ON.nodes.items():
            self.Graph.add_edge(self.chamber.singlestranded, on, kind = 'on_nucleation') #ADDPROPERTY
            self.Graph.add_edge(on, on.zipping[1], kind = 'zipping')
            for i, z in enumerate(on.zipping[2:], start=2):
                self.Graph.add_edge(z, on.zipping[i-1], kind = 'zipping')
            self.Graph.add_edge(on.zipping[-1], self.chamber.duplex, kind = 'zipping-end')

        L, R = self.chamber.split_offcores()
        for i, (left, right) in enumerate(zip(L, R)):
            self.Graph.add_edge(self.chamber.singlestranded, left, kind = 'off_nucleation')
            self.Graph.add_edge(self.chamber.singlestranded, right, kind = 'off_nucleation')
            if i > 0: 
                self.Graph.add_edge(left, L[i-1], kind = 'sliding')
                self.Graph.add_edge(right, R[i-1], kind = 'sliding')
        try: self.Graph.add_edge(L[-1], list(D.nodes())[0], kind = 'sliding-end')
        except IndexError: 
            if verbose: print('no left slidings as you can see from the empty list:', L)
            else: pass
        try: self.Graph.add_edge(R[-1], list(D.nodes())[0], kind = 'sliding-end')
        except IndexError: 
            if verbose: print('no right slidings as you can see from the empty list:', R)
            else: pass


    def clean_duplicates(self):
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



####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################



class Kinetics(object):
    def __init__(self, model: Model, kinetwork: Kinetwork, geometry: Geometry):
        
        self.model = model
        self.geometry = geometry
        self.T = model.kelvin
        self.space = model.space_dimensionality
        self.s1 = kinetwork.s1
        self.s2 = kinetwork.s2 
        
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
        self.georate = (np.power((self.geometry.theta/360),2))*(np.power((self.geometry.phi/360),2))*self.dlrate
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

    def set_zippingrate(self):
        """ first approximation is to take zipping
            equal to diffusion limited collision rate """
        self.zippingrate = self.dlrate



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
    def generalforward(self, kind):
        correspondence = {'sliding': self.slidingrate,
                          'sliding-end': self.slidingrate,
                          'zipping': self.zippingrate,
                          'zipping-end': self.zippingrate}
        return correspondence[kind]