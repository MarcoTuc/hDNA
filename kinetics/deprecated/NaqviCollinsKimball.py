from values import *
from kineticsfunctions import * 
import math
import numpy as np 
import scipy as sp 
from scipy import special, integrate


def front_term(D,R,beta):
    return phys['Na']*2*math.pi*D*R/(beta*100)

def J_den(u, D, beta):
    term1 = math.pi * u * special.j1(u) / (2 * np.sqrt(D/beta))
    term2 = 0.5 * math.pi * np.sqrt(D/beta) * special.j0(u)
    return (term1 + term2)**2

def Y_den(u, D, beta):
    term1 = math.pi * u * special.y1(u) / (2 * np.sqrt(D/beta))
    term2 = 0.5 * math.pi * np.sqrt(D/beta) * special.y0(u)
    return (term1 + term2)**2

def integrand(u, D, R, beta):
    alpha = D/(R**2)
    return (1/u)*np.exp(-alpha*(u**2))/((Y_den(u,D,beta)+J_den(u,D,beta)))

def naqvi_integral(D, R, beta):
    return integrate.quad(integrand, -0.0013, np.inf, args=(D,R,beta))

def naqvi(D,R,beta):
    return front_term(D,R,beta)*naqvi_integral(D,R,beta)[0]

