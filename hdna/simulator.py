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
from itertools import pairwise
from scipy.stats import expon

from .kinetwork import Kinetwork, Kinetics
from .model import Model, Options
from .strand import Strand



class Simulator(object):
    
    def __init__(self, model: Model, s1: Strand, s2: Strand, options=Options(), clean=False):
        
        self.model = model
        self.options = options
        self.kinet = Kinetwork(self.model, s1, s2)
        self.Graph = self.kinet.DG
        self.overview = self.kinet.overview

        methods = {"direct": jl.Direct()}
        self.method = methods[self.options.method]
        self.kinetics = self.kinet.kinetics
        self.biosim = jl.Network("biosim")

        self.sss = self.kinet.simplex #refers to simplex structure
        self.initialamount = self.options.initialamount

        if not clean:
            self.add_species()
            self.add_reactions()
            self.BSGraph()
        
            self.duplexindex = np.where(np.array([self.lt(str(e)) for e in list(self.biosim.species_list)]) == self.kinet.duplex)[0][0]
            self.duplex = self.kinet.duplex

            self.DIR = f'./{self.options.results_dir}/simulations/{self.options.stranditer}_{self.kinet.s1.sequence}'

    def add_species(self, verbose=False):
        """DONE"""
        # self.biosim <= jl.Species(self.tl(str(list(self.Graph.nodes())[0])), self.initialamount)
        ss = self.tl(self.sss)
        self.biosim <= jl.Species(ss, self.initialamount)
        for node in list(self.Graph.nodes())[1:]:
            if verbose: print(f'trying {node}')
            if verbose: print('TRANSLATION',self.tl(node))
            self.biosim <= jl.Species(self.tl(node))

    def add_reactions(self, verbose=False):
        for i, (e1, e2, data) in enumerate(self.Graph.edges.data()):
            state = self.Graph[e1][e2]['state']
            name = f"{state}-{i}"
            rate = data['k']
            etl1 = self.tl(e1)
            etl2 = self.tl(e2)
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


    def simulation(self):
        """ run the simulation and return the results """
        return jl.simulate(self.biosim, self.method, tfinal = self.options.runtime)
        # return {'c':jl.hcat(*sim.u).__array__(), 't':sim.t} OLD TRANSLATION TO DICTIONARY


    def ensemble(self):
        state, model = jl.parse_model(self.biosim)
        bar = tqdm(range(self.options.Nsim), leave = False)
        bar.set_description('Simulating')
        sim = [jl.simulate(state, model, self.method, tfinal = self.options.runtime) for _ in bar]
        if not self.options.make_sim_csv:
            return sim
        else:
            DIR_TRAJ = f'{self.DIR}/trajectories'
            try: os.makedirs(DIR_TRAJ)
            except FileExistsError: pass
            for run in sim:
                self.get_trajectory(run, weightlift=True, savetraj=False)
            for i, s in enumerate(sim[::int(len(sim)/self.options.trajstosave)]):
                traj = self.get_trajectory(s, weightlift=False, savetraj=True,) 
                traj.to_csv(f'{DIR_TRAJ}/run{i+1}.csv')
            self.options.stranditer += 1
            return sim
    

    def directsimulation(self):

        DIR_TRAJ = f'{self.DIR}/trajectories'
        DIR_TRAJ_FAIL = f'{self.DIR}/trajectories/failed'
        try:
            os.makedirs(DIR_TRAJ)
        except FileExistsError: pass
        try:
            os.makedirs(DIR_TRAJ_FAIL)
        except FileExistsError: pass
        fpts = []
        failed = 0
        state, model = jl.parse_model(self.biosim)
        bar = tqdm(range(self.options.Nsim), leave = False)
        bar.set_description('Simulating')

        for i in bar:
            sim = jl.simulate(state, model, self.method, tfinal = self.options.runtime)
            traj = self.get_trajectory(sim, weightlift=True, savetraj=True)
            if self.model.space_dimensionality == '3D':
                kcoll = self.kinetics.collisionrate
                collisiontime = expon(scale=1/kcoll).rvs()
            elif self.model.space_dimensionality == '2D':
                nuc = self.lt(traj[1])
                pos = self.kinet.DG.nodes[nuc]['obj'].sdist
                if pos <= self.kinetics.treshold:
                    kcollbot = self.kinetics.collisionbot
                    collisiontime  = expon(scale=1/kcollbot).rvs()
                else:
                    kcolltop = self.kinetics.collisiontop
                    collisiontime  = expon(scale=1/kcolltop).rvs()
            try:
                duplexationtime = sim.t[jl.findfirst(jl.isone, sim[self.duplexindex,:])-1]
                time = duplexationtime+collisiontime
                fpts.append(time)
                if i % int(self.options.Nsim/30) == 0:
                    pd.DataFrame(traj).to_csv(f'{DIR_TRAJ}/run{i+1}.csv')
            except TypeError:
                fpts.append(self.options.runtime)
                pd.DataFrame(traj).to_csv(f'{DIR_TRAJ_FAIL}/run{i+1}.csv')
                failed += 1

        print(f"{failed} simulations didn't produce a duplex.")
        print(f"That's {100*failed/self.options.Nsim}% of simulations")
        self.overview['failed'] = failed
        self.overview['fail%'] = 100*failed/self.options.Nsim
        self.options.stranditer += 1        
        return fpts, 1/np.mean(fpts)
    

    def get_trajectory(self, simulation, weightlift=True, savetraj=False): #TODO REDO IT MORE EFFICIENTLY (?)
        
        """ Create a method for appending a simulation trajectory to each simulation result
            Over than being pretty this is a good way to see inside simulations what's happening
            in order to see if the code is making stuff that makes actual physical sense.
            Checking this kind of stuff on the very big and sparse dataframes is unwise """

        if savetraj: self.trajectory = []
        states = list(self.Graph.nodes())
        # substitute 1 instead of 2 into the singlestranded states 
        for i, e in enumerate(simulation[0,:]): 
            if e == 2: simulation[0,:][i] = 1
        # Append the corresponding structure to the trajectory list
        if not weightlift:
            for step in simulation:
                index = jl.findall(jl.isone, step)[0] - 1
                if savetraj: self.trajectory.append(states[index])
        else:
            for i, (step, next) in enumerate(pairwise(simulation)):
                i_step = jl.findall(jl.isone, step)[0] - 1
                i_next = jl.findall(jl.isone, next)[0] - 1
                if savetraj: self.trajectory.append(states[i_step])
                try: self.BSG[states[i_step]][states[i_next]]['weight'] += 1
                except KeyError: pass
                if states[i_next] == self.duplex:
                    if savetraj: self.trajectory.append(states[i_next])
                    break
        # else:
        #     for i, (step, next) in enumerate(pairwise(simulation)):
        #         if i == 0:
        #             i_step = jl.findall(jl.isone, step)[0] - 1
        #             i_next = jl.findall(jl.isone, next)[0] - 1
        #         else:
        #             i_step = i_next
        #             i_next = jl.findall(jl.isone, next)[0] - 1
        #         if savetraj: self.trajectory.append(states[i_step])
        #         self.BSG[states[i_step]][states[i_next]]['weight'] += 1
        #         if states[i_next] == self.duplex:
        #             if savetraj: self.trajectory.append(states[i_next])
        #             break
            # Check for void trajectories - DISABLED 
            # if len(index) != 1: 
            #     if len(index) == 0: 
            #         # Cannot have void states (strand always exists)
            #         raise TrajectoryError(f'void state detected at step {step}')
            #     else:     
            #         # Cannot have more than one state at the same time 
            #         raise TrajectoryError(f'simultaneous states detected in step {step} at indices {[*list(index)]}')
            # index = index[0] - 1
        # if not self.options.rates_info:
        return self.trajectory
        # FRE = []
        # rates = []
        # names = []
        # for i, step in enumerate(self.trajectory[:-1], start=1):
        #     FRE.append(self.BSG.nodes[str(step)]['obj'].G)
        #     rates.append('{:e}'.format(self.BSG.edges[str(step),str(self.trajectory[i])]['rate']))
        #     names.append(self.BSG.edges[str(step),str(self.trajectory[i])]['name'])
        # DF = pd.DataFrame([self.trajectory, FRE, rates, names], 
        #                     index=['trajectory', 'DG', 'next step rate', 'step name']).T
        # return DF

    def BSGraph(self, verbose=False):
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
        self.BSG = nx.from_pandas_edgelist(dataframe, source='reactants', target='products', 
                                            edge_attr=['name', 'rate'], create_using=nx.DiGraph())
        # loop to inherit node properties from self.Graph to self.BSG
        for g in list(self.Graph.nodes):
            if verbose: print(g)
            for key in list(self.Graph.nodes[g].keys()):
                if verbose: print(self.Graph.nodes[g][key])
                self.BSG.nodes[g][key] = self.Graph.nodes[g][key]
        # loop to remove reaction number and get name as the original reaction state 
        # and add an empty weight 
        for edge in self.BSG.edges.data():
            edge[2]['state'] = edge[2]['name'].split('-')[0]
            edge[2]['weight'] = 1

        return self.BSG
    
    def save_graph(self, path=None, both=False):
        #convert node object to string of object type
        BSGsave = self.BSG.copy()
        if both: DGsave = self.kinet.DG.copy()
        for n in BSGsave.nodes.data():
            del n[1]['obj']
            try: del n[1]['tdx'] 
            except KeyError: pass
        if both:
            for n in DGsave.nodes.data():
                del n[1]['obj']
                try: del n[1]['tdx'] 
                except KeyError: pass
        if path != None:
            self.graphsaveformat(BSGsave)
            nx.write_gexf(BSGsave,f'{path}/{self.options.stranditer}_{self.kinet.s1.sequence}_graph_S.gexf')
            if both: 
                self.graphsaveformat(DGsave)
                nx.write_gexf(DGsave,f'{path}/{self.options.stranditer}_{self.kinet.s1.sequence}_graph_K.gexf')
        else:
            if self.options.graphsalone == 'own_folder': PATH = f'{self.options.results_dir}/graphs'
            elif self.options.graphsalone == 'strand_folder': PATH = self.DIR
            try: os.makedirs(PATH)
            except FileExistsError: pass 
            self.graphsaveformat(BSGsave)
            nx.write_gexf(BSGsave,f'{PATH}/{self.options.stranditer}_{self.kinet.s1.sequence}_graph_S.gexf')
            if both: 
                self.graphsaveformat(DGsave)
                nx.write_gexf(DGsave,f'{PATH}/{self.options.stranditer}_{self.kinet.s1.sequence}_graph_K.gexf')
            # except:
            #     self.BSGraph() 
            #     for n in self.BSG.nodes.data():
            #         n[1]['obj'] = str(type(n[1]['obj']))
            #     try: os.makedirs(PATH)
            #     except FileExistsError: pass 
            #     self.graphsaveformat(self.BSG)
            #     nx.write_gexf(self.BSG,f'{PATH}/{self.options.stranditer}_{self.kinet.s1.sequence}_graph_S.gexf')
            #     if both: nx.write_gexf(self.kinet.DG,f'{path}/{self.options.stranditer}_{self.kinet.s1.sequence}_graph_K.gexf')
        
                
    def mfpt(self, ensemble):
        fpts = self.fpts(ensemble)
        return np.mean(fpts)

    def fpts(self, ensemble):
        fpts = []
        failed = []
        for sim in ensemble:
            index = jl.findfirst(jl.isone, sim[self.duplexindex,:])
            try: fpts.append(sim.t[index-1])
            except TypeError: 
                fpts.append(self.options.runtime)
                failed.append(index)
        print(f"{len(failed)} simulations didn't produce a duplex.")
        print(f"That's {100*len(failed)/len(ensemble)}% of simulations")
        self.overview['failed'] = len(failed)
        self.overview['fail%'] = 100*len(failed)/len(ensemble)
        return fpts

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

    def graphsaveformat(self, Graph):
        for s, data in Graph.nodes.data():
            for key in data.keys():
                try: Graph.nodes[s][key] = '{:.3f}'.format(data[key])
                except ValueError: pass

        for e1, e2, data in Graph.edges.data():
            for key in data.keys():
                try: Graph[e1][e2][key] = '{:.3e}'.format(float(data[key]))
                except ValueError: pass

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


class TrajectoryError(Exception):
    def __init__(self, message):
        super().__init__(message)
