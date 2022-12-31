from values import *
from kineticsfunctions import * 
import math
import numpy as np 
import scipy as sp 
from scipy import special, integrate


def front_term(D):
    return phys['Na']*8*D/(math.pi*100)

def integrand(u, D, R):
    alpha = D/(R**2)
    return np.exp(-alpha*(u**2))/(u*(special.y0(u)**2 + special.j0(u)**2))

def naqvi_integral(D, R):
    return integrate.quad(integrand, 0, np.inf, args=(D,R))

def naqvi(D,R):
    return front_term(D)*naqvi_integral(D,R)[0]

