
import math
import numpy as np

import sys
sys.path.insert(0, '..')
sys.path.insert(0, '.')
from lidakinetics.values import *


#### THERMODYNAMICS 

def kelvin(temperature):
    return temperature + 273.15

def gibbs_free(enthalpy, entropy, temperature): #[DG]=Kcal/mol
    entropy_kcal = entropy/1000 #[S]=cal/mol*K
    return (enthalpy - (temperature)*entropy_kcal) # [H]=Kcal/mol 
    
def k_equilibrium(DG, temperature):
    return math.exp(-(DG/(phys['R(kcal/molK)']*(temperature))))   # (Kcal/mol)/((Kcal/mol*K)*K) -> adimensional




#### EYRING

def eyring_backwards(DG,Z,temperature):
    return (phys['k_boltz']*kelvin(temperature)/phys['h_planck'])*Z*math.exp(-(DG/(phys['R(kcal/molK)']*(kelvin(temperature)))))



#### 3D Diffusion Einstein Smoluchowski Relationships

def einsmol_spherical(radius, viscosity, temperature, units='cm'):
    if units == 'cm':
        return (phys['k_boltz']*temperature)/(6*math.pi*viscosity*radius)*10000
    if units == 'm':
        return (phys['k_boltz']*temperature)/(6*math.pi*viscosity*radius)

def smolu_diffusionlimited_equalparticles(viscosity, temperature):
    return (1/1000)*phys['Na']*(8/3)*phys['k_boltz_cm']*temperature/viscosity

def smolu_diffusionlimited_general(DA,DB,RA,RB):
    return 4*math.pi*(DA+DB)*(RA+RB)*phys['Na']/1000




#### Chew2019 rate function for two-dimensional activation limited reactions

#TODO
def Cc_chew2019_SCK(D, R):
    pass

def chew2019_al_MLM(time, diffusion, probability, voxel):
    """
    Activation-limited rate from chew2019
    This still has some strange stuff since 
    the rate doesn't change when changing probability
    """
    gamma = 0.57722
    b1 = 2*np.pi/np.sqrt(3)
    beta = b1*((1/probability) - 1)
    z = beta + np.log(48*time*diffusion/np.power(voxel,2))
    h2 = (np.power(gamma,2) - np.power(np.pi,2)/6)/np.power(z,3)
    Kchew = 4*np.pi*diffusion*((1/z) - gamma*(1/np.power(z,2)) - h2)
    return Kchew

def noyes_twodimensional(time, diffusion, probability, voxel):
    """
    Diffusion-limited rate from Noyes Theory
    as explained in Thorney-McConnell 1983
    """
    gamma = 0.57722
    beta = (1/probability)*np.pi*(1-probability)
    z = beta + np.log(16*diffusion*time*(1/np.power(voxel,2)))
    higher1 = gamma*(1/np.power(z,2))
    higher2 = ((1/6)*np.power(np.pi,2) - np.power(gamma,2))*(1/np.power(z,3))
    Kn = 4*np.pi*diffusion*((1/z) - higher1 - higher2) 
    return Kn



#### Geometric Rate 

def DNA_geometric_rate(DL_rate, azimutal_angle, longitudinal_angle):
    return (np.power((azimutal_angle/360),2))*(np.power((longitudinal_angle/360),2))*DL_rate


#### Collision Theory

# def collision_rate_3d(cross,mass):
#     return cross*np.power((4*phys['']))





#### DEPRECATED two dimensional diffusion-limited rates

# def tmc_2d_diffusionlimited(D,R,time):
#     return 4*math.pi*D*(math.log(4*D*time/(R**2)) - 2*phys['gamma'])*phys['Na']*(1/100)

# def melo5_diffusionlimited(D,R,beta):
#     return phys['Na']*2*math.pi*D*R/(beta*100)