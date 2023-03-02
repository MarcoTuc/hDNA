import nupack as nu 

class Model(object):

    def __init__(   self, 
                    sliding=None,
                    zipping=None,
                    sliding_filter=None,
                    material='dna', 
                    space_dimensionality='3D',
                    vescicleradius=0.5e-5, #decimeters (around a bacteria)
                    min_nucleation=1, 
                    sliding_cutoff=100, 
                    normalizeback=False,
                    stacking='nostacking', 
                    Na=1.0, 
                    Mg=0.0, 
                    celsius=26,
                    standard=False):
        """ I will use this class for passing experimental conditions """

        if standard:
            self.setgeometry(120, 270)
            self.setparams(zipping=4e7, sliding=1e7, sliding_filter=3)
       
        if material not in ['dna', 'rna']:
            raise ValueError("Material must be set to 'dna' or 'rna'")
        else: self.material = material 

        if space_dimensionality not in ['3D', '2D']:
            raise ValueError("Supported dimensionalities: 2D and 3D")
        else: self.space_dimensionality = space_dimensionality
        
        self.vescicleradius = vescicleradius

        allowedstackings = ['stacking', 'dangle-stacking', 'coaxial-stacking', 'nostacking']
        if stacking not in allowedstackings:
            raise ValueError(f'nupack only supports these stackings:\n{allowedstackings}')
        else: self.stacking = stacking
  
        self.min_nucleation = min_nucleation
        self.sliding_cutoff = sliding_cutoff
        self.normalizeback  = normalizeback
        
        self.celsius = celsius
        self.kelvin = celsius + 273.15
       
        self.Na = Na
        self.Mg = Mg 

        self.sliding = sliding
        self.zipping = zipping

        #do not touch 
        self.alpha = 1
        self.gamma = 0
        self.kappa = 1

        self.sliding_filter = sliding_filter

        self.nupack = nu.Model(material=self.material, 
                                    ensemble=self.stacking, 
                                    celsius=self.celsius, 
                                    sodium=self.Na, 
                                    magnesium= self.Mg) 

        if standard:
            self.setgeometry(120, 270)
            self.setparams(zipping=4e7, sliding=1e7, sliding_filter=3)

    def setparams(self, **kwargs):
        for i in kwargs.keys():
            if i not in ['zipping', 'sliding', 'sliding_filter']:
                raise ValueError('argument not valid')
        self.__dict__.update(kwargs)

    def setgeometry(self, theta=120, phi=270):
        self.theta = theta
        self.phi = phi


class Options(object):
    def __init__(self, 

                method="direct", 
                runtime=4e-6, 
                Nsim=1000,
                make_sim_csv=True,
                rates_info = True,
                save_graph_html = True,
                trajstosave=40,
                results_dir = './results', #beware if results_dir exists 
                graphsalone = 'strand_folder', 
                stranditer = 1
                ):
        
        # JULIA SIMULATION OPTIONS
        self.initialamount = 2 #never change this
        self.runtime = runtime
        self.Nsim = Nsim
        methods = ["direct"]
        if method in methods:
            self.method = method
        else: 
            message = f'method {method} is not yet implemented in the simulator. \n Choose from {methods}'
            raise NotImplementedError(message)

        # DATASAVING OPTIONS
        self.make_sim_csv = make_sim_csv
        self.rates_info = rates_info
        self.save_graph_html = save_graph_html
        self.results_dir = results_dir 
        self.graphsalone = graphsalone
        self.trajstosave = trajstosave
        self.stranditer = stranditer

        self.simtqdmdisable = True


# DEPRECATED
# class Geometry(object):

#     """ Here all the physical values needed for 
#     calculations will be defined and called for
#     other classes to be used """

#     def __init__(self, azimutal_angle, longitudinal_angle):
#         self.theta  = azimutal_angle
#         self.phi    = longitudinal_angle