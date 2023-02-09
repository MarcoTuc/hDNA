import csv 
import os

import networkx as nx
import pandas as pd 
import numpy as np

import juliacall
# print("Julia has",juliacall.Main.seval("Threads.nthreads()"),"threads")

jl = juliacall.newmodule("hDNA")

jl.seval('using BioSimulator')

from tqdm import tqdm 
from .kinetwork import Kinetwork, Kinetics
from .model import Model, Geometry

class Options(object):
    def __init__(self, 

                method="direct", 
                runtime=1e-5, 
                Nsim=500,

                make_sim_csv=True,
                rates_info = True,
                save_graph_html = True,
                trajstosave=10,
                results_dir = './results', #beware if results_dir exists 
                graphsalone = 'own_folder', 
                stranditer = 1
                
                ):
        
        # JULIA SIMULATION OPTIONS
        self.initialamount = 2
        self.runtime = runtime
        self.Nsim = Nsim
        methods = {"direct": jl.Direct()}
        self.method = methods[method]

        # DATASAVING OPTIONS
        self.make_sim_csv = make_sim_csv
        self.rates_info = rates_info
        self.save_graph_html = save_graph_html
        self.results_dir = results_dir 
        self.graphsalone = graphsalone
        self.trajstosave = trajstosave
        self.stranditer = stranditer


 #add some error here if the method is not a julia object or whatever 

class Simulator(object):
    
    def __init__(self, model: Model, kinetwork: Kinetwork, options=Options()):
        
        self.model = model 
        self.kinet = kinetwork
        self.Graph = kinetwork.Graph
        self.overview = self.kinet.overview

        self.kinetics = kinetwork.kinetics

        self.options = options 
        self.biosim = jl.Network("biosim")

        self.sss = self.kinet.chamber.singlestranded.structure
        self.initialamount = self.options.initialamount
        self.add_species()
        self.add_reactions()


    def add_species(self, verbose=False):
        """DONE"""
        # self.biosim <= jl.Species(self.tl(str(list(self.Graph.nodes())[0])), self.initialamount)
        ss = self.tl(self.kinet.chamber.singlestranded.structure) 
        self.biosim <= jl.Species(ss, self.initialamount)
        for node in list(self.Graph.nodes())[1:]:
            if verbose: print(f'trying {node}')
            if verbose: print('TRANSLATION',self.tl(node)) 
            self.biosim <= jl.Species(self.tl(node))

    def add_reactions(self, verbose=False):
        #NUCLEATIONS:
        for i, (e1, e2, data) in enumerate(self.Graph.edges.data()):
            name = f"{data['kind']}_{i}"
            rate = data['k']
            etl1 = self.tl(e1)
            etl2 = self.tl(e2)
            if verbose:
                print(name)
                print(rate)               
            if e1 == self.sss:
                rule = f"{etl1} + {etl1} --> {etl2}"
                if verbose: print(rule)
                self.biosim <= jl.Reaction(str(name), float(rate), str(rule))
            elif e2 == self.sss:
                rule = f"{etl1} --> {etl2} + {etl2}"
                if verbose: print(rule)
                self.biosim <= jl.Reaction(str(name), float(rate), str(rule))
            else: 
                rule = f"{etl1} --> {etl2}"
                if verbose: print(rule)
                self.biosim <= jl.Reaction(str(name), float(rate), str(rule))   
                

    
    def add_reactions_old(self):
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
        print(neighbors)
        sliding_factor = 1/(self.kinet.s1.length + self.kinet.s2.length - 1)
        zipping_factor = 1/(self.kinet.overview['on_nucleation'])
        for n in neighbors:
            i += 1
            """ JULIA PATHOLOGY 
                Can't use dots as names for reactions because
                Julia has a pathology with dots and parentheses 
                the most reliable notation is: 
                |   . --> o     |
                |   ( --> b     |
                |   ) --> d     |
                |   + --> A     |                   """
            neigh = self.tl(n)
            state = self.Graph.nodes[n]['state'] 
            factor = sliding_factor if state == 'off_register' else zipping_factor
            # FORWARD BIMOLECULAR NUCLEATION
            name = f"f_nucleation{i}"; rule = f"{ss} + {ss} --> {neigh}"
            kf = self.kinetics.georate * factor
            self.biosim <= jl.Reaction(name, kf, rule)
            self.Graph.edges[SS,n]['kf'] = kf   # Also add the rate as a new property to the graph, it will be useful for visualization

            # BACKWARD UNIMOLECULAR DISSOCIATION            
            name = f"b_nucleation{i}"; rule = f"{neigh} --> {ss} + {ss}"
            DG = self.Graph.nodes[n]['obj'].G
            kb = self.kinetics.k_back(kf, DG)
            self.biosim <= jl.Reaction(name, kb, rule)
            self.Graph.edges[SS,n]['kb'] = kb   # Also add the rate as a new property to the graph, it will be useful for visualization

        # now make a graph without the singlestranded state for all subsequent transitions 
        subgraph = self.Graph.copy()
        subgraph.remove_node(SS)
        
        for n1, n2, data in list(subgraph.edges.data()):
            i += 1
            
            # if data['kind'] in ['sliding', 'sliding-end']:
            l1 = self.Graph.nodes[n1]['obj'].total_nucleations
            l2 = self.Graph.nodes[n2]['obj'].total_nucleations
            if l1 < l2:             # Left slidings
                n1o = self.tl(n1)
                n2o = self.tl(n2)
            else:                   # Right slidings
                n1o = self.tl(n2)
                n2o = self.tl(n1)
            # else:
            #     n1o = self.tl(n1)
            #     n2o = self.tl(n2)

            # FORWARD SLIDINGS AND ZIPPINGS
            name = f"f_{data['kind']}_{i}"; rule = f"{n1o} --> {n2o}"
            kf = data['k']
            self.biosim <= jl.Reaction(name, kf, rule)
            self.Graph.edges[n1,n2][f'kf'] = kf   # Also add the rate as a new property to the graph, it will be useful for visualization

            # BACKWARD SLIDINGS AND ZIPPINGS 
            name = f"b_{data['kind']}_{i}"; rule = f'{n2o} --> {n1o}'
            DG = self.Graph.nodes[n2]['obj'].G - self.Graph.nodes[n1]['obj'].G
            if self.Graph.nodes[n1]['obj'].total_nucleations > self.Graph.nodes[n2]['obj'].total_nucleations: DG = - DG
            kb = data['k']
            self.biosim <= jl.Reaction(name, kb, rule)
            self.Graph.edges[n1,n2][f'kb'] = kb   # Also add the rate as a new property to the graph, it will be useful for visualization

            """ Note that for slidings coming from the right like: .......((((+.......)))) <--> ........(((+........))) 
                the forward rate will always denote the passage from the strand with less basepairs with the one with more.
                By now the network visualization is kinda broken but who cares, rates are where they should be. 
                This happens because of how nodes are ordered inside kinetwork's Graph """

    def simulation(self):
        """ run the simulation and return the results """
        return jl.simulate(self.biosim, self.options.method, tfinal = self.options.runtime)
        # return {'c':jl.hcat(*sim.u).__array__(), 't':sim.t} OLD TRANSLATION TO DICTIONARY 
    
    def ensemble(self):
        state, model = jl.parse_model(self.biosim)
        sim = [jl.simulate(state, model, self.options.method, tfinal = self.options.runtime) for _ in tqdm(range(self.options.Nsim))]
        if not self.options.make_sim_csv:
            return sim
        else:
            DIR = f'./{self.options.results_dir}/simulations/{self.options.stranditer}_{self.kinet.s1.sequence}'
            DIR_TRAJ = f'{DIR}/trajectories'  
            try: os.makedirs(DIR_TRAJ)
            except FileExistsError: pass 
            for i, s in enumerate(sim[::int(len(sim)/self.options.trajstosave)]):
                traj = self.get_trajectory(s) 
                traj.to_csv(f'{DIR_TRAJ}/run{i+1}.csv')
            if self.options.graphsalone == 'own_folder': self.save_graph(f'{self.options.results_dir}/graphs')
            elif self.options.graphsalone == 'strand_folder': self.save_graph(f'{DIR}')
            self.options.stranditer += 1
            return sim


    def mfpts(self, ensemble):
        """
        node referring to Duplex should be the last one so the result spot for duplex should just be something like simresults[-1] i hope 
        first i need to understand how the simresult looks like out of the juliacalled biosimlator """

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
        self.overview['failed'] = len(failed)
        self.overview['fail%'] = 100*len(failed)/len(ensemble)
        return np.mean(fpts)
    

    def get_trajectory(self, simulation): #TODO REDO IT MORE EFFICIENTLY (?)
        
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
        if not self.options.rates_info:
            return self.trajectory
        else:
            G = self.DiGraph()
            rates = []
            names = []
            for i, step in enumerate(self.trajectory[:-1], start=1):
                rates.append('{:e}'.format(G.edges[str(step),str(self.trajectory[i])]['rate']))
                names.append(G.edges[str(step),str(self.trajectory[i])]['name'])
            DF = pd.DataFrame([self.trajectory, rates, names], 
                                index=['trajectory', 'next step rate', 'step name']).T
            return DF

    def DiGraph(self, verbose=False):
        """ Return a digraph post-biosim to see if it matches with the pre-biosim graph.
            This is useful to see if there are any errors in the kinetwork-biosim translation"""
        R = self.biosim.reaction_list
        properties = ['name', 'rate', 'reactants', 'products']
        names = []; rates = []; reags = []; prods = []
        for r in R:
            names.append(str(r))
            rates.append(R[r].rate)
            reags.append(self.lt(str(list((R[r].reactants))[0])))
            prods.append(self.lt(str(list((R[r].products))[0])))
        dataframe = pd.DataFrame([names, rates, reags, prods], index=properties).T
        self.digraph = nx.from_pandas_edgelist(dataframe, source='reactants', target='products', 
                                            edge_attr=['name', 'rate'], create_using=nx.DiGraph())
        # loop for inheriting node properties from self.Graph to self.digraph
        for g in list(self.Graph.nodes):
            if verbose: print(g)
            for key in list(self.Graph.nodes[g].keys()):
                # print(self.digraph.nodes[g][key])
                if verbose: print(self.Graph.nodes[g][key])
                self.digraph.nodes[g][key] = self.Graph.nodes[g][key]

        return self.digraph
    
    def save_graph(self, PATH):
        #convert node object to string of object type
        try: 
            for n in self.digraph.nodes.data():
                n[1]['obj'] = str(type(n[1]['obj']))
            try: os.makedirs(PATH)
            except FileExistsError: pass 
            nx.write_gexf(self.digraph,f'{PATH}/{self.options.stranditer}_{self.kinet.s1.sequence}_graph.gexf')
        except:
            self.DiGraph() 
            for n in self.digraph.nodes.data():
                n[1]['obj'] = str(type(n[1]['obj']))
            try: os.makedirs(PATH)
            except FileExistsError: pass 
            nx.write_gexf(self.digraph,f'{PATH}/{self.options.stranditer}_{self.kinet.s1.sequence}_graph_S.gexf')
        

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

    #########################
    ### Reporting Methods ###
    #########################

    def print_properties(self, *args):
        """ Define a method for printing the properties of nodes/edges 
            in the graph in a nice looking way."""
        pass

    def save_overview(self, path, name):
        df = pd.DataFrame.from_dict([self.overview]).T
        df.rename(columns={np.int64(0):'values'}, inplace=True)
        df.index.rename('states', inplace=True)
        df.to_csv(f'{path}/{name}.csv')

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
    #     for sim in tqdm(range(self.options.Nsim)):
    #         data = ensemble[sim].u.__array__()
    #         time = ensemble[sim].t
    #         new = {'c':data,'t':time}
    #         out.append(new)
    #     return new
    #     # return [{'c':jl.hcat(*ensemble[i].u).__array__(), 't':ensemble[i].t} for i in tqdm(range(len(ensemble)))]




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