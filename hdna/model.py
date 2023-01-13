class Model(object):

    def __init__(   self, 
                    material, 
                    space_dimensionality,
                    Na=1.0, 
                    Mg=0.0, 
                    celsius=26):
        """ I will use this class for passing experimental conditions """
       
        if material != ('dna' or 'rna'):
            raise ValueError("Material must be set to 'dna' or 'rna'")
        else: self.material = material 

        # if space_dimensionality != '2D' or '3D':
        #     raise ValueError("Allowed space dimensionality: 2D and 3D")

        if space_dimensionality not in ['3D', '2D']:
            raise ValueError("Supported dimensionalities: 2D and 3D")

        self.space_dimensionality = space_dimensionality
        
        self.celsius = celsius
        self.kelvin = celsius + 273.15
       
        self.Na = Na
        self.Mg = Mg 

class Geometry(object):

    """ Here all the physical values needed for 
    calculations will be defined and called for
    other classes to be used """

    def __init__(self, azimutal_angle, longitudinal_angle):
        self.theta  = azimutal_angle
        self.phi    = longitudinal_angle