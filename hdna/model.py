import nupack as nu 

class Model(object):

    def __init__(   self, 
                    material, 
                    space_dimensionality,
                    min_nucleation=1, 
                    sliding_cutoff=1, 
                    stacking='nostacking', 
                    sliding=2e7,
                    zipping=2e9,
                    Na=1.0, 
                    Mg=0.0, 
                    celsius=26):
        """ I will use this class for passing experimental conditions """
       
        if material not in ['dna', 'rna']:
            raise ValueError("Material must be set to 'dna' or 'rna'")
        else: self.material = material 

        if space_dimensionality not in ['3D', '2D']:
            raise ValueError("Supported dimensionalities: 2D and 3D")
        else: self.space_dimensionality = space_dimensionality

        allowedstackings = ['stacking', 'dangle-stacking', 'coaxial-stacking', 'nostacking']
        if stacking not in allowedstackings:
            raise ValueError(f'nupack only supports these stackings:\n{allowedstackings}')
        else: self.stacking = stacking
  
        self.min_nucleation = min_nucleation
        self.sliding_cutoff = sliding_cutoff
        
        self.celsius = celsius
        self.kelvin = celsius + 273.15
       
        self.Na = Na
        self.Mg = Mg 

        self.sliding = sliding
        self.zipping = zipping

        self.nupack = nu.Model(material=self.material, 
                                    ensemble=self.stacking, 
                                    celsius=self.celsius, 
                                    sodium=self.Na, 
                                    magnesium= self.Mg) 


class Geometry(object):

    """ Here all the physical values needed for 
    calculations will be defined and called for
    other classes to be used """

    def __init__(self, azimutal_angle, longitudinal_angle):
        self.theta  = azimutal_angle
        self.phi    = longitudinal_angle