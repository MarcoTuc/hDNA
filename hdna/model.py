import nupack as nu 

class Model(object):

    def __init__(   self, 
                    sliding=None,
                    zipping=None,
                    sliding_filter=None,
                    material='dna', 
                    space_dimensionality='3D',
                    min_nucleation=1, 
                    sliding_cutoff=100, 
                    stacking='nostacking', 
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
        self.sliding_filter = sliding_filter
        
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

    def setparams(self, zipping, sliding, sfilter=4):
        self.zipping = zipping
        self.sliding = sliding
        self.sliding_filter = sfilter

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
        self.initialamount = 2
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