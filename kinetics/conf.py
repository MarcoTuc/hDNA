import sys
sys.path.insert(0, 'lidakinetics')
sys.path.insert(0, 'nupackrelated')
sys.path.insert(0, 'dnapoly')
sys.path.insert(0, '..')
sys.path.insert(0, '.')

from lidakinetics.values import *
from lidakinetics.kineticsfunctions import *
from dnapoly import dnageometry
#from NaqviCollinsKimball import naqvi_integral

# GENERAL PARAMETERS
TEMPERATURE = kelvin(26)
VISCOSITY_BULK = hydroparams['H20_viscosity']['value']
HOURS = 0.1
TIME = 60*60*HOURS

# DNA GEOMETRIC PARAMETERS
ssDNA_RADIUS = dnageometry.ssDNA_GEOM['cylinder_radius_cm']

# DIFFUSION PARAMETERS
DIFFUSION3D = einsmol_spherical(ssDNA_RADIUS, VISCOSITY_BULK, TEMPERATURE,units='cm')
DIFFUSION2D = Diffusion2D_exp['Filippov04'] 

# print('3D Diffusion:',DIFFUSION3D)
# print('2D Diffusion:',DIFFUSION2D)
# print('{:3e}'.format(smolu_diffusionlimited_equalparticles(VISCOSITY_BULK,TEMPERATURE)))
# print('{:3e}'.format(smolu_diffusionlimited_general(DIFFUSION3D,DIFFUSION3D,ssDNA_RADIUS,ssDNA_RADIUS)))
# print('{:3e}'.format(tmc_2d_diffusionlimited(2*DIFFUSION2D,2*ssDNA_RADIUS,TIME)))
# print('{:3e}'.format(melo5_diffusionlimited(2*DIFFUSION2D,2*ssDNA_RADIUS,0.5)))


def sciprint(string:str):
    print('{:3e}'.format(string))
