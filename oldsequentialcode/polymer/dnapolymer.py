##############################################################
### Here I will put all the functions related to computing ###
### DNA related polymer physics properties                 ###
##############################################################
import sys

import numpy as np 
import lidakinetics.values as v
from dnapoly.dnageometry import *

def ssDNA_weight(sequence):
    nucleotide_weights = {'A':331.2218, 'T':322.2085, 'C':307.1971, 'G':347.2212}
    return sum([nucleotide_weights[sequence[i]] for i in range(len(sequence))])

def FJC_endtoend(n_monomers, mono_length, R_max_out=False):
    l2 = np.power(mono_length, 2)
    R_max = n_monomers*mono_length
    if R_max_out == False:
        return l2*n_monomers
    else:
        return l2*n_monomers, R_max

def FJC_endtoend_FloryDistribution(r,n,l):
    beta    = np.sqrt(3/2)/(np.sqrt(n)*l)
    coeff   = np.power((beta/np.sqrt(np.pi)),2)
    expo    = np.exp(-np.power(beta,2)*np.power(r,2))
    vol     = 4*np.pi*np.power(r,2) 
    return coeff*expo

def kuhn_transform(n,l,theta=0,C=7):
    R_max = n*l*np.cos(theta/2)
    b = (C*n*np.power(l,2))/R_max
    N = np.power(R_max,2)/(C*n*np.power(l,2))
    return N, b 

def get_kuhnchain_persistence(persistence, mono_length, n_monomers):
    b = 2*persistence
    N = (mono_length*n_monomers)/b
    return b, N

def wormlike_endtoend(lp, n_monomers, mono_length, R_max_out=False):
    '''
    lp:       persistence length
    n_monomers:     number of monomers in the chain 
    mono_length:    length of the monomers
    '''
    R_max = n_monomers*mono_length

    if R_max_out == False: 
        return 2*lp*R_max - 2*(lp**2)*(1 - np.exp(-(R_max/lp)))
    else: 
        return 2*lp*R_max - 2*(lp**2)*(1 - np.exp(-(R_max/lp))), R_max


def zimm_model_D(N_kuhn, b_kuhn, viscosity, temperature):
    prefactor = 8/(3*np.sqrt(6*np.power(np.pi,3)))
    thermal_energy = (v.phys['k_boltz_cm']*(temperature))
    radius_zimm = b_kuhn * np.sqrt(N_kuhn) * 1e-7
    return prefactor * thermal_energy / (viscosity * radius_zimm)


def px_realchain(x):
    return 0.278*(x**0.28)*np.exp(-1.206*(x**2.43))
def px_idealchain(x):
    return ((3/2*np.pi)**(3/2))*np.exp(-1.5*(x**2))


def gyration_radius_FJC(n_monomers):
    lp = 2.223
    ld = 0.676
    lambda_0 = (lp*ld)/3 
    return np.sqrt(n_monomers*lambda_0)