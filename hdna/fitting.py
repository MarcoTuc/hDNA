import pandas as pd
import numpy as np 
from scipy.stats import gamma, expon

def gammafit(fpts):
    alpha, loc, beta = gamma.fit(fpts)
    return gamma(a=alpha, loc=loc, scale=beta)

def expfit(fpts):
    loc, scale = expon.fit(fpts)
    return expon(loc=loc, scale=scale)