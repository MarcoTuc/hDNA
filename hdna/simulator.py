import csv 
import pickle 

import networkx as nx
import pandas as pd 
import numpy as np

import juliacall
# print("Julia has",juliacall.Main.seval("Threads.nthreads()"),"threads")

jl = juliacall.Main

jl.seval('using BioSimulator')

from tqdm import tqdm 
from hdna.kinetwork import Kinetwork, Kinetics
from hdna.model import Model, Geometry

class Options(object):
    def __init__(self, method="direct", runtime=1e-5, Nmonte=500):
        self.runtime = runtime
        self.Nmonte = Nmonte

        methods = {"direct": jl.Direct()}

        self.method = methods[method] #add some error here if the method is not a julia object or whatever 

class Simulator(object):
    
    def __init__(self, model: Model, kinetwork: Kinetwork, kinetics: Kinetics, initialamount=2, options=Options()):
        
        self.model = model 
        self.kinet = kinetwork
        self.Graph = kinetwork.Graph

        self.kinetics = kinetics 

        self.options = options 
        self.biosim = jl.Network("biosim")

        self.initialamount = initialamount
        self.add_species()
        self.add_reactions()


    def add_species(self):
        """DONE"""
        # self.biosim <= jl.Species(self.tl(str(list(self.Graph.nodes())[0])), self.initialamount)
        ss = self.tl(self.kinet.chamber.singlestranded.structure) 
        self.biosim <= jl.Species(ss, self.initialamount)
        for node in list(self.Graph.nodes())[1:]:
            self.biosim <= jl.Species(self.tl(node))


    def add_reactions(self):
        i = 0  #for now reaction names will just be sequential numbers
        """PSEUDOCODE:
        
        for node in self.Graph:
            get neighborhood of node 
            pairs = node.pairs (how many base pairs the node has)
            state = node.state (what is it)
            for node in neighborhood:
                self.biosim <= jl.Reaction(reaction_name, reaction_rate, reaction_string (use the f"{stuff}" method of Giovanni))
        """
        # start from the neighborhood of the singlestranded state to add nucleations 
        SS = self.kinet.chamber.singlestranded.structure 
        ss = self.tl(SS)
        neighbors = list(nx.neighbors(self.Graph, SS))
        sliding_factor = 1/len(neighbors)
        for n in neighbors:
            i += 1
            """ JULIA PATHOLOGY 
                Can't use dots as names for reactions because
                Julia has a pathology with dots and also 
                parenthesis (fuck) the most reliable notation is: 
                |   . --> o     |
                |   ( --> b     |
                |   ) --> d     |
                |   + --> A     |                   """
            neigh  = self.tl(n)

            # FORWARD BIMOLECULAR NUCLEATION
            name = f"f_nucleation{i}"; rule = f"{ss} + {ss} --> {neigh}"
            kf = self.kinetics.georate * sliding_factor
            self.biosim <= jl.Reaction(name, kf, rule)
            self.Graph.edges[SS,n]['kf'] = kf   # Also add the rate as a new property to the graph, it will be useful for visualization

            # BACKWARD UNIMOLECULAR DISSOCIATION            
            name = f"b_nucleation{i}"; rule = f"{neigh} --> {ss} + {ss}"
            DG = self.Graph.nodes[n]['object'].G
            kb = self.kinetics.k_back(kf, DG)
            self.biosim <= jl.Reaction(name, kb, rule)
            self.Graph.edges[SS,n]['kb'] = kb   # Also add the rate as a new property to the graph, it will be useful for visualization

        # now make a graph without the singlestranded state for all subsequent transitions 
        subgraph = self.Graph.copy()
        subgraph.remove_node(SS)
        
        for n1, n2, data in list(subgraph.edges.data()):
            i += 1                       # TODO RATES (data will be used to get rates etc)
            
            if data['kind'] in ['sliding', 'sliding-end']:
                l1 = self.Graph.nodes[n1]['object'].total_nucleations
                l2 = self.Graph.nodes[n2]['object'].total_nucleations
                if l1 < l2:             # Left slidings
                    n1o = self.tl(n1)
                    n2o = self.tl(n2)
                else:                   # Right slidings
                    n1o = self.tl(n2)
                    n2o = self.tl(n1)
            else:
                n1o = self.tl(n1)
                n2o = self.tl(n2)

            # FORWARD SLIDINGS AND ZIPPINGS
            name = f"f_{data['kind']}_{i}"; rule = f"{n1o} --> {n2o}"
            kf = self.kinetics.generalforward(data['kind'])
            self.biosim <= jl.Reaction(name, kf, rule)
            self.Graph.edges[n1,n2][f'kf'] = kf   # Also add the rate as a new property to the graph, it will be useful for visualization

            # BACKWARD SLIDINGS AND ZIPPINGS 
            name = f"b_{data['kind']}_{i}"; rule = f'{n2o} --> {n1o}'
            DG = self.Graph.nodes[n2]['object'].G - self.Graph.nodes[n1]['object'].G
            if self.Graph.nodes[n1]['object'].total_nucleations > self.Graph.nodes[n2]['object'].total_nucleations: DG = - DG
            kb = self.kinetics.k_back(kf, DG)
            self.biosim <= jl.Reaction(name, kb, rule)
            self.Graph.edges[n1,n2][f'kb'] = kb   # Also add the rate as a new property to the graph, it will be useful for visualization

            """ Note that for slidings coming from the right like: .......((((+.......)))) <--> ........(((+........))) 
                the forward rate will always denote the passage from the strand with less basepairs with the one with more.
                By now the network visualization is kinda broken but who cares, rates are where they should be. 
                This happens because of how nodes are ordered inside kinetwork's Graph """

    def simulation(self):
        """ run the simulation and return the results """
        return jl.simulate(self.biosim, self.options.method, tfinal = self.options.runtime)
        # return {'c':jl.hcat(*sim.u).__array__(), 't':sim.t}
    
    def ensemble(self):
        state, model = jl.parse_model(self.biosim)
        return [jl.simulate(state, model, self.options.method, tfinal = self.options.runtime) for _ in tqdm(range(self.options.Nmonte))]

    def mfpts(self, ensemble):
        """
        node referring to Duplex should be the last one so the result spot for duplex should just be something like simresults[-1] i hope 
        first i need to understand how the simresult looks like out of the juliacalled biosimlator 

        PSEUDOCODE:

        def condition(x):
            return 1 == y (when the value of duplex gets to one for the firs time we get the mfpt)
        def findmfpt(sim):
            index_r = findfirst(condition, result[:,duplexindex])
            firstpassage = etc etc 

        mfpts = []
        time = simresults.time
        for sim in simresults:
            mfpts.append(self.findmfpt(sim)) make another method that just gets the mfpt for a single simulation result 
        """

        fpts = []
        failed = []
        for sim in ensemble:
            index = jl.findfirst(jl.isone, sim[-1,:])
            try: fpts.append(sim.t[index-1])
            except TypeError: 
                fpts.append(self.options.runtime)
                failed.append(index)
        print(f"{len(failed)} simulations didn't produce a duplex.")
        print(f"That's {100*len(failed)/len(ensemble)}% of simulations")
        return np.mean(fpts)
        
        
        ############ DEPRECATED MFPTS ROUTINES ################
        # fpts = []
        # failed = []
        # try: data = self.ensemble_to_dict(ensemble)
        # except: data = ensemble 
        # for i, sim in enumerate(data):
        #     fpi = np.argmax(sim['c'][-1] == 1)
        #     fpts.append(data['t'][fpi])
        # print(f"{len(failed)} simulations didn't produce a duplex.")
        # print(f"That's {100*len(failed)/len(data)}% of simulations")
        # return np.mean(fpts)
    
        # fpts = []
        # failed = []
        # try: data = self.ensemble_to_dict(ensemble)
        # except: data = ensemble 
        # duplex = list(self.Graph.nodes())[-1]
        # for i, sim in enumerate(data):
        #     try: 
        #         index = np.where(sim[duplex] == 1)[0][0]
        #         fpt = sim['time'][index]
        #         fpts.append(fpt)
        #     except: 
        #         # print(f"simulation at index {i} has not duplexed")
        #         failed.append(i)
        #         fpts.append(self.options.runtime)
        #         # print(f"appended runtime as fpt for failed simulation {i}: {self.options.runtime}")
        # # print(f"list of failed simulations: {failed}")
        # print(f"{len(failed)} simulations didn't produce a duplex.")
        # print(f"That's {100*len(failed)/len(data)}% of simulations")
        # return np.mean(fpts)
    

    def get_trajectory(self, simulation): #TODO REDO IT EFFICIENTLY
        
        """ Create a method for appending a simulation trajectory to each simulation result
            Over than being pretty this is a good way to see inside simulations what's happening
            in order to see if the code is making stuff that makes actual physical sense.
            Checking this kind of stuff on the very big and sparse dataframes is unwise """
        
        self.trajectory = []
        states = list(self.Graph.nodes())
        # substitute 1 instead of 2 into the singlestranded states 
        for i, e in enumerate(simulation[0,:]): 
            if e == 2: simulation[0,:][i] = 1
        # Append the corresponding structure to the trajectory list
        for step in range(len(simulation)-1):
            stepvector = simulation[:,step]
            index = jl.findall(jl.isone, stepvector)
            # Check for errors
            if len(index) != 1: 
                if len(index) == 0: 
                    # Cannot have void states (strand always exists)
                    raise TrajectoryError(f'void state detected at step {step}')
                else:     
                    # Cannot have more than one state at the same time 
                    raise TrajectoryError(f'simultaneous states detected in step {step} at indices {[*list(index)]}')
            index = index[0] - 1
            self.trajectory.append(states[index])
        
        return self.trajectory




        ###### DUE TO CHANGES IN ENSEMBLE TO DICT IT DOESN'T WORK ANYMORE AT THE MOMENT
        # trajectory = []
        # try: data = simulation.drop('time', axis=1)
        # except KeyError: data = simulation
        # for index, row in data.iterrows():
        #     trajectory.append(row[row!=0].index[0])
        # return trajectory
        

    ######################
    ### Helper Methods ###
    ######################

    def tl(self, s):
        """ Used to fix julia pathology """
        table = str.maketrans({'.':'o', '(':'b', ')':'d', '+':'A'})
        return s.translate(table)
    
    def lt(self,s):
        table = str.maketrans({'o':'.', 'b':'(', 'd':')', 'A':'+'})
        return s.translate(table)

    ########################
    ### Printing Methods ###
    ########################

    def print_properties(self, *args):
        """ Define a method for printing the properties of nodes/edges 
            in the graph in a nice looking way."""
        pass

    ########################
    ### Property Methods ###
    ########################

    @property
    def duplex(self):
        return list(self.Graph.nodes())[-1]



class TrajectoryError(Exception):
    def __init__(self, message):
        super().__init__(message)







# ############################ DEPRECATED SHIT #######################################

#     def julia_csv_pandas(self, ensemble):
#         out = []
#         for sim in ensemble:
#             with open('translator.csv', 'w+') as file:
#                 file.truncate(0)
#                 writer = csv.writer(file)
#                 writer.writerows(sim)
#                 df = pd.read_csv('translator.csv', names=list(self.Graph.nodes()))
#             out.append(df)
#         return out 


#     def julia_pickle_pandas(self, ensemble):
#         out = []
#         for sim in ensemble:
#             file = open('mbare', 'wb')
#             pickle.dump(sim, file)
#             file.close()
#             file = open('mbare', 'rb')
#             df = pd.read_pickle(file)
#             file.close()
#             out.append(df)
#         return out 



    # def ensemble_to_dataframe(self, ensemble):
    #     """ convert the ensemble of simulations to a pandas dataframe """
    #     out = []
    #     for sim in ensemble:
    #         df = pd.DataFrame(sim.u, columns=list(self.Graph.nodes()))
    #         df['time'] = sim.t
    #         out.append(df)
    #     return out

    # def ensemble_to_dict(self, ensemble):
    #     out = []
    #     for sim in tqdm(range(self.options.Nmonte)):
    #         data = ensemble[sim].u.__array__()
    #         time = ensemble[sim].t
    #         new = {'c':data,'t':time}
    #         out.append(new)
    #     return new
    #     # return [{'c':jl.hcat(*ensemble[i].u).__array__(), 't':ensemble[i].t} for i in tqdm(range(len(ensemble)))]
