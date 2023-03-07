import numpy as np 

from typing import NamedTuple

from .model import Model
from .strand import Strand
from .params import *



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

        if self.model.space_dimensionality == '3D':
            ### Default compute relevant rates
            self.diffusionlimited()
            self.geometric_rate()
        elif self.model.space_dimensionality == '2D':
            self.georate = self.surfgeorate()


    """ ------------- ZIPPING RELATED -----------------"""

    def set_zippingrate(self, setrate):
        """ first approximation is to take zipping
            equal to diffusion limited collision rate """
        if setrate == 'diffusionlimited': self.zippingrate = self.dlrate
        else: self.zippingrate = setrate

    
    """ ------------- SLIDING RELATED -----------------"""

    def set_slidingrate(self, sliding): #NOTE NEED TO EXPAND THIS TO TRY ALL RANGES
        self.slidingrate = sliding
    
    def unif_scaling(self, nucleations): # --> TODO
        return 1/nucleations

    def z_scaling(self, nucleations): # --> TODO
        """
        Returns:
        - the nucleation partition function as sum over nodes of e^dg_node/KbT
        - list of boltzmann weights per node (use a dictionary)
        Then boltzmann weighting will just be indexing the given node
        from the dictionary and dividing the value by Z 
        """
        pass

    def gammasliding(self, dgs):
        return self.model.alpha * np.exp( self.model.gamma + (self.model.kappa * ((dgs) / (self.phys['R(kcal/molK)'] * (self.T)))))
        # 1 / ( 1 + np.exp( self.model.gamma + (dgs / (self.phys['R(kcal/molK)'] * (self.T))))) #HERTELGAMMASLIDING

    def overlappingspheres(self, d, r1, r2):
        separation = r1+r2-d
        if separation <= 0:
            return 0
        else:
            term1 = separation**2
            term2 = (d**2 + 2*d*r1 - 3*(r1**2) + 2*d*r2 + 6*r1*r2 - 3*(r2**2))
            overlap = np.pi*term1*term2/(12*d)
            sphere1  = (4/3)*np.pi*(r1**3)
            sphere2  = (4/3)*np.pi*(r2**3)
            return overlap/(sphere1+sphere2)

    """ ------------- KINETICS-THERMODYNAMICS RELATION METHODS -----------------"""
    
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
    

    """ ------------- THREE-DIMENSIONAL DIFFUSION LIMITED NUCLEATIONS -----------------"""

    #3333333333333333333333333333333333333333333333333333333333333333333333333333333
    # >>>>>>>>------------------- 3D DIFFUSION -----------------------<<<<<<<<<<<<<<

    def einsmol_spherical(self, radius):
        """ returns einstein smoluchowski diffusivity for 
            spherical particles in (cm^2)/s units """
        return self.phys['k_boltz_cm']*self.T/(6*np.pi*self.viscosity*radius)

    ########################## Vanilla
    def ss_strands_size(self):
        self.size1 = self.s1.length*self.simplex_geometry['basepairdistance']
        self.size2 = self.s2.length*self.simplex_geometry['basepairdistance']
    def vanilla_diffusivities(self):
        self.ss_strands_size()
        self.vD1 = self.einsmol_spherical(self.size1)
        self.vD2 = self.einsmol_spherical(self.size2)
    
    ##################### Polymer Physics
    def gyradiuses(self):
        lambda_0 = (self.simplex_geometry['persistence_length']*self.simplex_geometry['basepairdistance'])/3
        self.gr1 = np.sqrt(self.s1.length*lambda_0)
        self.gr2 = np.sqrt(self.s2.length*lambda_0)
        return self.gr1, self.gr2

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
    
    # Check dnapolymer.py inside ./oldsequentialcode/polymer 
    # for all the other polymer physics related functions

    # >>>>>>>>------------------- 3D DIFFUSION LIMITED RATES -----------------------<<<<<<<<<<<<<<
    def diffusionlimited(self, kind='ppi'):
        """ Smoluchowski 1916 classical formula """
        if kind == 'ppi':  
            self.ppi_diffusivities()
            cm2dmcubic = 1e-3 #dimensional correction from cubic cm to cubic dm to get 1/(M*s) kinetic rates
            self.dlrate = self.phys['Na']*4*np.pi*(self.pD1+self.pD2)*(self.gr1+self.gr2) * cm2dmcubic
            return self.dlrate
        elif kind == 'vanilla':
            self.vanilla_diffusivities()
            dc = 1e-3 #dimensional correction from cubic cm to cubic dm to get 1/(M*s) kinetic rates
            self.dlrate = self.phys['Na']*4*np.pi*(self.vD1+self.vD2)*(self.size1+self.size2) * dc
            return self.dlrate
        else: raise ValueError(f'{kind} not implemented')

    def bulksteric(self):
        return (np.power((self.model.theta/360),2))*(np.power((self.model.phi/360),2))

    def geometric_rate(self):
        self.georate = self.bulksteric() * self.dlrate
        return self.georate 

    def closedconfscaling(self, p_circular): # --> TODO 
        """ p_circular => probability of circularization 
            circularization will impede nucleation because ... 
            (just make a picture of it in your mind for now) """ 
        self.georate = self.georate * (1 - p_circular)
        #TODO implement methods to calculate p_circular from 
        #end-to-end probability distributions from polymer physics 


    """ ------------- TWO DIMENSIONAL DIFFUSION LIMITED NUCLEATIONS -----------------"""

    #2222222222222222222222222222222222222222222222222222222222222222222222222222222
    # >>>>>>>>------------------- 2D DIFFUSION -----------------------<<<<<<<<<<<<<<
    
    def tdiff_saffdelb(self):
        """ Saffman-Delbruck Membrane Translational Diffusion"""
        e = MMGEO.LIP_RAD/HYDRO.LSD
        return ((CONST.KB*self.model.kelvin)/(4*np.pi*HYDRO.ETA_MEM)) * (np.log(2/e) - CONST.GAMMA)
    
    def surfacesteric(self):
        return np.power((self.model.theta/360),2)

    def surfgeorate(self):
        return self.ksphere_sano() * self.surfacesteric()

    # >>>>>>>>------------------- 2D DIFFUSION LIMITED KINETICD -----------------------<<<<<<<<<<<<<<
    
    def ksphere_sano(self):
        # right now just approximate 2D gyradius with 3D one
        a, _ = self.gyradiuses()
        b    = 5e-5 #took 1 micrometer as the radius, idk 
        surf = 4*np.pi*np.power(b*1e-1,2)
        r = (a + b)/a 
        D = self.tdiff_saffdelb()
        coeff1 = (np.power(a + b, 2)/D)
        coeff2 = 2/(1-(np.power(1/r, 2)))
        coeff3 = np.log(r)
        return coeff1 * ( coeff2 * coeff3 - 1 ) * surf * CONST.NA

    def kc_chew2019(self, time):
        """
        Activation-limited rate from chew2019
        This still has some strange stuff since 
        the rate doesn't change when changing probability
        """
        D = self.tdiff_saffdelb()
        b1 = 2*np.pi/np.sqrt(3)
        beta = b1*((1/self.surfacesteric()) - 1)
        z = beta + np.log(48*time*D/np.power(sum(self.gyradiuses()),2))
        h2 = (np.power(CONST.GAMMA,2) - np.power(np.pi,2)/6)/np.power(z,3)
        Kchew = 4*np.pi*D*((1/z) - CONST.GAMMA*(1/np.power(z,2)) - h2)
        return Kchew * 1e-2 * CONST.NA
    def kc_noyes(self, time):
        """
        Diffusion-limited rate from Noyes Theory
        as explained in Thorney-McConnell 1983
        """
        D = self.tdiff_saffdelb()
        beta = (1/self.surfacesteric())*np.pi*(1-self.surfacesteric())
        z = beta + np.log(16*D*time*(1/np.power(sum(self.gyradiuses()),2)))
        higher1 = CONST.GAMMA*(1/np.power(z,2))
        higher2 = ((1/6)*np.power(np.pi,2) - np.power(CONST.GAMMA,2))*(1/np.power(z,3))
        Kn = 4*np.pi*D*((1/z) - higher1 - higher2) 
        return Kn * 1e-2 * CONST.NA
    def kc_torney(self, time):
        D = self.tdiff_saffdelb()
        ln = np.log(4*D*time/np.power(sum(self.gyradiuses()),2))
        return 4*np.pi*D*(1/(ln - 2*CONST.GAMMA)) * 1e-2 * CONST.NA

    """ ------------- EQUILIBRIUM AND BACKWARD RATES -----------------"""

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
    
    
    def avgunzip(self):
        return self.zippingrate*np.exp((-2)/(CONST.R*self.model.kelvin))

