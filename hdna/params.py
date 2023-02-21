from typing import NamedTuple
from dataclasses import dataclass
import PhysicalQuantities as pq 
from PhysicalQuantities import q

@dataclass(frozen=True)
class Constants:
    KB      = 1.380649e-19 # J/K = N*m/K = m^2*kg/s^2*K = 10^4 cm^2*kg/s^2*K
    NA      = 6.023e23     # mol^-1
    R       = 1.987e-3     # Kcal/mol*K
    GAMMA   = .57722       # 

@dataclass(frozen=True)
class Hydroparams:
    MU_H2O  = 0.8701e-5         # kg/cm*s
    ETA_MEM = 5e-10             # Pa*s*m = kg/s
    LSD = ETA_MEM/(2*MU_H2O)    # Saffman-Delbruck Length - centimeters

@dataclass(frozen=True)
class DuplexGeo: #centimeters
    MONODIST    = 0.34e-7
    PERSISTENCE = 100*MONODIST
    CY_RADIUS   = 2e-7

@dataclass(frozen=True)
class SimplexGeo: #centimeters
    MONODIST    = 0.676e-7
    PERSISTENCE = 2.223e-7
    CY_RADIUS   = 1e-7

@dataclass(frozen=True)
class MembraneGeo: #centimeters
    LIP_RAD     = 2e-7 

HYDRO = Hydroparams()
CONST = Constants()
DXGEO = DuplexGeo()
SXGEO = SimplexGeo()
MMGEO = MembraneGeo()

