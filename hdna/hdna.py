import os 
import sys 
from tqdm import tqdm

import pandas as pd 
import numpy as np 

from scipy.optimize import basinhopping

from hdna import *

OPTDIR = 'results'

class HDNA:
    def __init__(self, data, expname, model=Model(), options = Options()):
        
        self.data = data
        self.model = model
        self.options = options
        self.iteration = 1
        sequences = self.data['sequences']
        expvalues = self.data['experimental']
        self.tupledata = [(s, e) for s, e in zip(sequences, expvalues)]
        self.RUN = 1
        self.EXPNAME = expname 
        self.folder = f'./{OPTDIR}/{self.EXPNAME}'
        self.check_dir()
    
    def set_options(self, *args):
        self.options = Options(*args)


    def run(self, rates):
        
        savedata = self.data
        
        self.model.setparams(zipping=rates[0], sliding=rates[1])
        self.model.setgeometry()
        savefolder = f'{self.folder}/run_{self.RUN}'
        self.make_dir(savefolder)

        saveparams = {'zipping': rates[0],
                      'sliding': rates[1]}
        pd.DataFrame(saveparams, index=[0]).T.to_csv(f"{savefolder}/parameters.csv")

        f = open(f'{savefolder}/console.txt', 'w')
        sys.stdout = Tee(f)

        i = 0
        bar = tqdm(self.tupledata, leave=False)
        for seq, exp in bar:
            bar.set_description(f'{seq}')
            i += 1
            print(f'Strand number {i}: {seq}')
            A = Strand(self.model, seq)
            self.options.results_dir = savefolder
            self.options.stranditer = i
            print('Embedding network into biosimulator network model...')
            S = Simulator(self.model, A, A.complementary(), self.options)
            for d in S.overview.items():
                print(f'{d[0]} = {d[1]}')
            fptimes, modrate = S.directsimulation()
            print(f'{len(fptimes)} fpts have been produced')
            pd.DataFrame(fptimes).to_csv(f"{S.DIR}/fpts.csv")
            savedata.loc[seq, 'computational'] = modrate
            for col, val in zip(S.overview.keys(), S.overview.values()):
                savedata.loc[seq, col] = val
            self.save_plots(fptimes, S, seq, exp, modrate)
            S.save_graph()
            print(f"experimental rate: {'{:e}'.format(float(exp))}")
            print(f"computed rate:     {'{:e}'.format(float(modrate))}", '\n')


        #Routine for LMSE output
        for i, x in savedata[['experimental', 'computational']].iterrows():
            savedata.loc[i, 'error'] = np.power(float(x['experimental']) - float(x['computational']), 2)  
        savedata.to_csv(f"{savefolder}/simulationdata.csv")
        valplot(savedata, f'scatter_{self.EXPNAME}_run-{self.RUN}', writepath=savefolder, theme='dark')
        self.RUN += 1

        return savedata['error'].sum()

    
    def optimize(self):
        initial = [2e7, 2e7]
        optimizer = basinhopping(self.run, initial, niter=5, stepsize=1e7, T=300)

    
    def save_plots(self, fpts, simulatore, sequence, exp, mod):
            try: 
                fitgamma = gammafit(fpts)
                histotime(fpts, fitgamma, simulatore.options.runtime, exp=exp, mod=mod, name=f'{sequence}gamma_histo', writepath=simulatore.DIR, theme='dark')
            except: pass
            try: 
                fitexp = expfit(fpts)
                histotime(fpts, fitexp, simulatore.options.runtime, exp=exp, mod=mod, name=f'{sequence}_exp_histo', writepath=simulatore.DIR, theme='dark') 
            except: pass
            percomplot(fpts, writepath=simulatore.DIR, name=f'{sequence}_percentplot')

    def make_dir(self, RESULTS_DIR):
         # Directory Check 
        if os.path.isdir(RESULTS_DIR): 
            i = 0
            while True: 
                i += 1
                permission = input(f'Folder {RESULTS_DIR} already exists, do you want to overwrite old experiments? [Y,N] ')
                if permission.lower().startswith('y'):
                    print('>>>> overwriting old simulations')
                    break
                elif permission.lower().startswith('n') or i == 3:
                    asknew = input(f'Do you want to make a new folder? [Y,N] ')
                    if asknew.lower().startswith('y'):
                        control = False
                        while True: 
                            NEWEXPNAME = input(f'Write the new experiment name to put inside ./results folder:\n')
                            RESULTS_DIR = f"results/{NEWEXPNAME}"
                            if NEWEXPNAME != self.EXPNAME: 
                                os.makedirs(RESULTS_DIR)
                                self.EXPNAME = NEWEXPNAME
                                break
                            else: print('Bro, put a new name not the old one... \n')
                        break 
                    elif permission.lower().startswith('n'): 
                        print(">>>> stopping the program")
                        sys.exit()
                print("yes or not?") 
        else:
            os.makedirs(RESULTS_DIR)
    
    def check_dir(self):
        self.make_dir(self.folder)





    # def run(self, rates):
        
    #     savedata = self.data
        
    #     self.model.setparams(zipping=rates[0], sliding=rates[1])
    #     self.model.setgeometry()
    #     savefolder = f'{self.folder}/run_{self.RUN}'
    #     self.make_dir(savefolder)

    #     saveparams = {'zipping': rates[0],
    #                   'sliding': rates[1]}
    #     pd.DataFrame(saveparams, index=[0]).T.to_csv(f"{savefolder}/parameters.csv")

    #     f = open(f'{savefolder}/console.txt', 'w')
    #     sys.stdout = Tee(f)

    #     i = 0
    #     bar = tqdm(self.tupledata, leave=False)
    #     for seq, exp in bar:
    #         bar.set_description(f'{seq}')
    #         i += 1
    #         print(f'Strand number {i}: {seq}')
    #         A = Strand(self.model, seq)
    #         B = A.complementary()
    #         self.options.results_dir = savefolder
    #         self.options.stranditer = i
    #         print('Embedding network into biosimulator network model...')
    #         S = Simulator(self.model, A, B, self.options)
    #         for d in S.overview.items():
    #             print(f'{d[0]} = {d[1]}')
    #         trajectories = S.ensemble()
    #         fptimes = S.fpts(trajectories)
    #         modrate = 1/(np.mean(fptimes))
    #         savedata.loc[seq, 'computational'] = modrate
    #         for col, val in zip(S.overview.keys(), S.overview.values()):
    #             savedata.loc[seq, col] = val
    #         self.save_plots(fptimes, S, seq, exp, modrate)
    #         print(f"experimental rate: {'{:e}'.format(float(exp))}")
    #         print(f"computed rate:     {'{:e}'.format(float(modrate))}", '\n')


    #     #Routine for LMSE output
    #     for i, x in savedata[['experimental', 'computational']].iterrows():
    #         savedata.loc[i, 'error'] = np.power(float(x['experimental']) - float(x['computational']), 2)  
    #     savedata.to_csv(f"{savefolder}/simulationdata.csv")
    #     valplot(savedata, f'scatter_{self.EXPNAME}_run-{self.RUN}', writepath=savefolder, theme='dark')
    #     self.RUN += 1

    #     return savedata['error'].sum()