import nupack as nu 
import numpy as np 
import pandas as pd 
import networkx as nx 
import re

from itertools import pairwise, tee

from .strand import Strand
from .model import Model
from .kinetwork import Kinetics


class Complex(object):
    
    def __init__(   self, 
                    model: Model, 
                    s1: Strand, 
                    s2: Strand, 
                    state,
                    structure=None,
                    dpxdist=None,
                ):
        
        if type(model) != Model:
            raise TypeError("Model must be an instance of hdna.Model")
        self.model = model

        if type(s1) != Strand:
            raise TypeError("s1 must be an instance of hdna.Strand")
        if type(s2) != Strand:
            raise TypeError("s2 must be an instance of hdna.Strand")

        self.s1 = s1    #53
        self.s2 = s2    #35

        self.l1 = s1.length
        self.l2 = s2.length

        # Think of this as self.initialstructure
        # I didn't change it because of laziness 
        self.structure = structure

        self.dpxdist = dpxdist

        self.zipgraph = nx.DiGraph()

        self.getnupackmodel()
        self.getnupackproperties()
        self.G = self.structureG(self.structure)

        self.mismatches = []
        self._get_mismatches()

        self.possible_states = [ None,
                            'singlestranded',
                            'duplex',
                            'zipping',
                            'on_nucleation',
                            'off_nucleation', 
                            'backfray',
                            'sliding' ]   

    
        if state not in self.possible_states:
            raise ValueError(f'state must be one among {self.possible_states} but you gave {state}')
        else: 
            self.state = state 
        
        if self.state == 'duplex':
            #TODO update this when considering mismatches 
            self.consecutive_nucleations = min(self.l1, self.l2)
            self.total_nucleations = self.consecutive_nucleations
            self.structure = '('*self.s1.length+'+'+')'*self.s2.length
            self.dpxdist = 0
        elif self.state == 'singlestranded':
            self.total_nucleations = 0
            self.consecutive_nucleations = 0
            self.structure = '.'*self.s1.length+'+'+'.'*self.s2.length
            self.G = 0
        else:
            self.total_nucleations = self.totbasepairs(self.structure)
            self.consecutive_nucleations = self.maxconsbp(self.structure)

        """ These variables are used to say if the current
        Coils object is duplexed, is offregister or is onregister
        to be comunicated to classes higher in the hierarchy"""

        if self.state == 'on_nucleation': self.zipping_trajectory()
    
    def set_state(self, target):
        if target != self.state:
            self.state = target

        
    
######################################
##### Self-calculated Properties #####
######################################

    #TODO: ADD ABASIC DEFECTS (not supported by nupack) 
    # ab = {
    # 'GUG':  {'A':-4.3,'C':-3,'G':-5.7,'T':-4.5},
    # 'CUC':  {'A':-7.6,'C':-11.3,'G':-8.3,'T':-11.3},
    # }

    def totbasepairs(self, structure):
        nl = structure.count('(')
        nr = structure.count(')')
        if nl == nr: tot = nl
        else: raise BrokenPipeError(f'left and right nucleation should match: {structure}')
        return tot
    
    def maxconsbp(self, structure):
        l, r = structure.split('+')
        maxl = len(max(l.split('.')))
        maxr = len(max(r.split('.')))
        if maxl == maxr: cons = maxl
        else: raise BrokenPipeError(f'left and right nucleation should match: {structure}')
        return cons


    def zippingraph(self):

        def get_i(lst):
            indices = []
            for i, (el1, el2) in enumerate(pairwise(lst)):
                if el1 == '.' and el2 != '.':
                    indices.append(i)
                if el1 != '.' and el2 == '.':
                    indices.append(i+1)
            return indices

        def update_structure(string, side, character: str):
            updated = string
            moveslist = []
            moves = get_i(string)[::side]
            for m in moves:
                up = updated[:m] + character + updated[m+1:]
                moveslist.append(up)
            return moveslist
    
        def leafs(left, right, graph):
            lmoves = update_structure(left, 1, '(')
            rmoves = update_structure(right, -1, ')')
            for lmove, rmove in zip(lmoves, rmoves):
                
                #FORWARD
                self.zipgraph.add_edge(Zippo(self.model, self.s1, self.s2, structure='+'.join([left, right])), 
                                       Zippo(self.model, self.s1, self.s2, structure='+'.join([lmove,rmove])))
                #BACKWARD
                self.zipgraph.add_edge(Zippo(self.model, self.s1, self.s2, structure='+'.join([lmove,rmove])),
                                       Zippo(self.model, self.s1, self.s2, structure='+'.join([left, right])))      
            
                yield list(leafs(lmove,rmove,graph))


    

    def zipping_trajectory(self):
        """
        CAUTION:    use only for native nucleation states with no double nucleation.
                    don't use with eg "((...((..+..))...)) strands
        UPDATE:     This has been corrected but still one 
                    should use caution against this method"""
        
        def get_i(lst):
            indices = []
            for i, el in enumerate(lst):
                if i > 0 and el != lst[i-1]:
                    indices.append(i)
            return indices

        def update_structure(string, character: str):
            indices = get_i(string)
            indices_inv = get_i(string[::-1])
            updated = string
            for index in indices:
                updated = updated[:index-1] + character + updated[index:]
            updated_inv = updated[::-1]
            for index in indices_inv:
                updated_inv = updated_inv[:index-1] + character + updated_inv[index:]
            return updated_inv[::-1]

        if self.state == 'on_nucleation': 
            self.zipping = []
            left, right = self.structure.split('+')
            while '.' in left and right:
                left = update_structure(left, '(')
                right = update_structure(right, ')')
                step = '+'.join([left,right])
                # The generated structure is parsed to see if it is correctly base-paired 
                candidatestep = self.parse_structure(step, self.s1, self.s2)
                self.zipping.append(Zippo(self.model, self.s1, self.s2, state='zipping', structure=candidatestep, dpxdist=0))
            self.zipping.pop(-1) #remove the spurious duplex element generated at the end of the while loop
            return self.zipping
        else: self.zipping = None; return self.zipping
        
    def inherit_zipping(self, listofzippings):
        self.zipping = listofzippings


#####################################
##### Nupack-dependant Methods  #####
#####################################

    def structureG(self, structure=None):
        """ Nupack related properties """
        if structure == None: structure = self.structure
        nuStrand1 = nu.Strand(self.s1.sequence, name = 'a')
        nuStrand2 = nu.Strand(self.s2.sequence, name = 'b') # An inversion here is needed because in this program strands are defined as 5-3 against 3-5 but in NUPACK all strands are defined 5-3 and the program takes care to turn them around and so on
        nuStructure = nu.Structure(structure)
        dG = nu.structure_energy(strands=[nuStrand1,nuStrand2], structure=nuStructure, model=self.nupackmodel)
        return dG

    def getnupackmodel(self):
        self.nupackmodel = nu.Model(material=self.model.material, 
                                    ensemble=self.model.stacking, 
                                    celsius=self.model.celsius, 
                                    sodium=self.model.Na, 
                                    magnesium= self.model.Mg) 

    def getnupackproperties(self):
            """ Nupack related properties """
            self.getnupackmodel()
            self.nuStrand1 = nu.Strand(self.s1.sequence, name = 'a')
            self.nuStrand2 = nu.Strand(self.s2.sequence, name = 'b') # An inversion here is needed because in this program strands are defined as 5-3 against 3-5 but in NUPACK all strands are defined 5-3 and the program takes care to turn them around and so on
            self.nuComplex = nu.Complex([self.nuStrand1,self.nuStrand2], name = 'c')
            j = nu.pfunc(self.nuComplex, self.nupackmodel)
            self.duplexG = j[1]
            self.duplexZ = float(j[0])
            self.mfe = nu.mfe([self.nuStrand1, self.nuStrand2], model=self.nupackmodel)


#############################################
##### General Structure Parsing Methods #####
#############################################
    
    """ This is a more general structure parsing routine 
        borrowed from the chamber class and needed here to
        parse zipping trajectories which for now don't account
        for wattson-crick base pairing.
        Since structures here have the normal dot-bracket notation
        this method doesn't work for ì structures. Maybe I will
        make a Structure class for handling all these structure
        related things """

    def parse_structure(self, structure, seq_a, seq_b):
        """
        Converts the internal general structures made with
        'ì' and '.' to dot-bracket structures considering the
        Wattson Crick pairings of bases
        """
        struct_a, struct_b = structure.split('+')
        seq_a = seq_a.sequence
        seq_b = seq_b.sequence 
        sx, dx = self.structurecut(seq_a, seq_b[::-1], struct_a, struct_b[::-1])
        patch_sx = ''.join(['\(' for i in range(len(sx))])
        patch_dx = ''.join(['\)' for i in range(len(dx))])
        out1 = re.sub(patch_sx,sx,struct_a)
        out2 = re.sub(patch_dx,dx,struct_b)
        structureout = out1+"+"+out2
        return structureout


    def structurecut(self, string1, string2, structure1, structure2):
        cut1 = ''.join([n1 for n1, s1 in zip(string1, structure1) if s1 == '('])
        cut2 = ''.join([n2 for n2, s2 in zip(string2, structure2) if s2 == ')'])
        sx = ''
        dx = ''
        # print(cut1, cut2)
        for n1, n2 in zip(cut1, cut2):
            if self._iswattsoncrick(n1, n2):
                sx += '('
                dx += ')'
            else: 
                sx += '.'
                dx += '.' 
        return sx, dx[::-1]

##########################
##### Helper Methods #####
##########################

    def _iswattsoncrick(self,n1,n2):
        """
        check for mismatches
        """
        if (n1 == 'A' and n2 == 'T') or (n1 == 'C' and n2 == 'G') or (n1 == 'T' and n2 == 'A') or (n1 == 'G' and n2 == 'C'):
            return True
        else:
            return False 

    def _get_mismatches(self):
        for nuc1, nuc2 in zip(self.s1.sequence, self.s2.sequence):
            self.mismatches.append(self._iswattsoncrick(nuc1, nuc2))

    def splitstructure(self):
        self.struct_a, self.struct_b = self.structure.split('+')
        return self.struct_a, self.struct_b

######################################
##### Methods Under Construction #####
######################################

    def correctedstrands(self):
        """ if given sequences have abasic defects this 
        method corrects the given Strands to non abasic 
        equivalent to pass them to nupack, the abasic
        penalty will be taken care of in the next method"""
        pass

    def abasicpenalty(self):
        """ Check the presence of an abasic defect
        and take the penalty out of the currently 
        calculated free energy """
        pass

######################
##### Properties #####
######################

    @property
    def sequences(self):
        return self.s1.sequence+"+"+self.s2.sequence



""" ---------------- SUBCLASSES ------------------ """



class Zippo(Complex):
    def __init__(self, 
                    model: Model, 
                    s1: Strand, 
                    s2: Strand, 
                    state,
                    structure,
                    dpxdist=None):
        super().__init__(model,s1,s2,state,structure)
        self.dpxdist = dpxdist


class Sliding(Complex):
    def __init__(self,
                    model:Model,
                    s1: Strand,
                    s2: Strand,
                    state,
                    structure,
                    dpxdist):
        super().__init__(model,s1,s2,state,structure)
        self.dpxdist = dpxdist #distance in terms of basepairs from the sliding to the duplex
        # self.backfraying_trajectory()
        self.off_nucleations()
        # self.sbackfray = [b.structure for b in self.backfray]
        # self.gbackfray = [b.G for b in self.backfray]
        # try: self.maxstable = self.backfray[np.argmin(self.gbackfray)]
        # except ValueError: self.maxstable = self 

    def off_nucleations(self, verbose=False):
        
        def get_ix(string, char):
            indices = []
            for i, e in enumerate(list(string)):
                if e == char:
                    indices.append(i)
            return indices 

        def trans(l):
            trans = str.maketrans({'(': 'b', ')': 'd'})
            return l.translate(trans)

        def backtrans(l):
            trans = str.maketrans({'b': '(', 'd': ')'})
            return l.translate(trans)
        
        def transparens(l):
            trans = str.maketrans({'(':'.', ')':'.'})
            return l.translate(trans)

        def nwise(iterable,n):
            iterators = tee(iterable, n)
            for i, iter in enumerate(iterators):
                for j in range(i):
                    next(iter, None)
            return zip(*iterators)

        def replace(index_list,character,string):
            string=list(string)
            for index in index_list:
                string[index]=character
            return "".join(string)

        def update(struct, verbose=False):
            ixl = get_ix(struct, '(')
            ixr = get_ix(struct, ')')
            ixb = get_ix(struct, 'b')
            ixd = get_ix(struct, 'd')
            if verbose: print(ixl, ixr, ixb, ixd)
            try: nldo = min([j for j in ixl if j<ixb[0]],  key=lambda x:abs(x-ixb[0]))
            except ValueError: nldo = ixb[0]
            try: nlup = min([j for j in ixl if j>ixb[-1]], key=lambda x:(abs(x-ixb[-1])))
            except ValueError: nlup = ixb[-1]
            try: nrdo = min([j for j in ixr if j<ixd[0]],  key=lambda x:abs(x-ixd[0]))
            except ValueError: nrdo = ixd[0]
            try: nrup = min([j for j in ixr if j>ixd[-1]], key=lambda x:abs(x-ixd[-1]))
            except ValueError: nrup = ixd[-1]
            if verbose: print(ixb);print(ixl);print(nldo, nlup, nrdo, nrup)
            struct = replace([nldo, nlup], 'b', struct)
            struct = replace([nrdo, nrup], 'd', struct)
            return struct 
        
        ixl = get_ix(self.structure, '(')
        ixr = get_ix(self.structure, ')')

        if self.consecutive_nucleations > self.model.min_nucleation:
            self.backfray = []
            for l, r in zip(
                nwise(ixl, self.model.min_nucleation),
                nwise(ixr[::-1], self.model.min_nucleation)):
                new = replace(l, 'b', self.structure)
                new = replace(r, 'd', new)
                newtrans = backtrans(transparens(new))
                offcore = Complex(
                    self.model, 
                    self.s1, 
                    self.s2, 
                    state='off_nucleation',
                    structure=newtrans,
                    dpxdist=self.dpxdist)
                self.backfray.append(offcore)
                itszippings = []
                while new != trans(self.structure):
                    new = update(new)
                    newtrans = backtrans(transparens(new))
                    if verbose: 
                        print('slidings backfray routine produced','\n',newtrans)
                        print('from source', self.structure)
                    itszippings.append(Zippo(
                        self.model, 
                        self.s1, 
                        self.s2, 
                        state='backfray',
                        structure=newtrans,
                        dpxdist=self.dpxdist))
                offcore.inherit_zipping(itszippings[:-1])
        else: self.backfray = []
            
        
        # #internal utility function for getting indices 
        # def get_ix(string, char):
        #     indices = []
        #     for i, e in enumerate(list(string)):
        #         if e == char:
        #             indices.append(i)
        #     return indices 
        # #internal utility function for string translation
        # def trans(l):
        #     trans = str.maketrans({'.': '.', '(': 'ì', ')': 'ì'})
        #     return l.translate(trans)
        # #internal utility function for generating nwise iterators 
        # #see itertools pairwise but now nwise dude 
        # def nwise(iterable,n):
        #     iterators = tee(iterable, n)
        #     for i, iter in enumerate(iterators):
        #         for j in range(i):
        #             next(iter, None)
        #     return zip(*iterators)
        # #iuf for replacement of string with character at list of indices 
        # def replace(index_list,character,string):
        #     string=list(string)
        #     for index in index_list:
        #         string[index]=character
        #     return "".join(string)
        
        # #Actual output of the method 
        # ixl = get_ix(self.structure, '(') #get left strand indices
        # ixr = get_ix(self.structure, ')') #get right strand indices 
        # self.offcores = []
        # for l, r in zip(
        #     nwise(ixl, self.model.min_nucleation), 
        #     nwise(ixr, self.model.min_nucleation)):
        #     self.offcores.append([*l, *r], 'ì', self.structure)

    def backfraying_trajectory(self):
        from itertools import zip_longest as zipp

        def get_ix(string, char):
            indices = []
            for i, e in enumerate(list(string)):
                if e == char:
                    indices.append(i)
            return indices 

        def backfray_structure():
            l, r = self.structure.split('+')
            lix = get_ix(l, '('); ll = int(len(lix)/2)
            rix = get_ix(r, ')'); rl = int(len(rix)/2)
            backfray = []
            for li, lil, ri, ril in zipp(
                lix[:ll], lix[ll:][::-1], 
                rix[:rl], rix[rl:][::-1]):
                # print(li, ri, lil, ril )
                try: l = l[:li] + '.' + l[li+1:]
                except TypeError: pass
                try: l = l[:lil] + '.' + l[lil+1:]
                except TypeError: pass
                try: r = r[:ri] + '.' + r[ri+1:]
                except TypeError: pass
                try: r = r[:ril] + '.' + r[ril+1:]
                except TypeError: pass 
                s = '+'.join([l,r])
                backfray.append(
                    Zippo(self.model, self.s1, self.s2, state='backfray', structure=s))
            backfray.pop(-1)
            return backfray
        
        if self.total_nucleations >= self.model.min_nucleation and self.consecutive_nucleations >= self.model.min_nucleation:
            steps = backfray_structure()
            self.backfray = [self]
            for step in steps:
                self.backfray.append(step)
            self.backfray = self.backfray[::-1]
            self.backfray[0].set_state('off_nucleation')
        else: 
            self.set_state('off_nucleation')
            self.backfray = []
    
    @property
    def bfempty(self):
        return self.backfray == []


        # if False: #part of backfraystructure
        #     print('mbarestructure',self.structure)
        #     indices = get_ix(string, char)
        #     # indices_inv = get_ix(string[::-1], character)
        #     print('indices dir',indices)
        #     # print('indices inv',indices_inv)
        #     updated = string
        #     for index in [indices[0], indices[-1]]:
        #         # print(string, index)
        #         print(f'updated {char}', index, updated)
        #         updated = updated[:index] + '.' + updated[index+1:]
        #         print(f'updated {char}', index, updated)
        #     # updated_inv = updated[::-1]
        #     # for index in [indices_inv[0], indices_inv[-1]]:
        #         # print(string, index)
        #         # updated_inv = updated_inv[:index-1] + '.' + updated_inv[index:]
        #     return updated

        # if False:
        #     if self.consecutive_nucleations > 3:
        #         backfray = []
        #         l, r = self.structure.split('+')
        #         while '(' in l and ')' in r: 
        #             l = backfray_structure(l, '(')
        #             r = backfray_structure(r, ')')
        #             step = '+'.join([l,r])
        #             # print(l)
        #             # print(r)
        #             # print(self.structure)
        #             # print(step)
        #             backfray.append(Zippo(self.model, self.s1, self.s2, step, self.structureG(step)))
        #     else: print("i'm a sweet boy:", self.structure)


################################################################################################################################
##########################################################À CEMETERY ###########################################################
############################################### MAY THESE FUNCTIONS REST IN PEACE ##############################################
################################################################################################################################ 
  
    # def gotnucleationhere(self, mincore):
    #     self.nucleation = None
    #     if self.l1 != self.l2:
    #         return ValueError('This method cannot be used for complexes made with different length strands')
    #     if self.nucleation_size >= mincore:
    #         self.nucleation == True
    #     else:
    #         self.nucleation = False

    # def maxnucleation(self):
    #     self.nucleation_size = 0
    #     for n in range(0, len(self.mismatches) + 1):
    #         if self.checknucleation(n):
    #             self.nucleation_size = n
    #     return self.nucleation_size

    # def checknucleation(self, mincore):
    #     self.nucleation = None
    #     if self.l1 != self.l2:
    #         return ValueError('This method cannot be used for complexes made with different length strands')
    #     for b in range(len(self.mismatches) - (mincore-1)):
    #         if all(self.mismatches[b:b+mincore]):
    #             return True 
    #     return False 

    # def nucleationsize(self):
    #     """ This method gives back two quantities:
    #         - Total number of non-mismatched nucleated base pairs
    #         - Size of the biggest consecutive sequence of nucleated base pairs """
    #     A, B = self.splitstructure()
    #     A, B = pd.Series(A), pd.Series(B)
    #     nA = A.str.count("\(").sum()
    #     nB = B.str.count("\)").sum()
    #     if nA != nB: raise ValueError('total number of based pairs should match for each strand')
    #     else: self.total_nucleations = nA 
    #     splitsA = A.str.split('[^(]')
    #     splitsB = B.str.split('[^)]')
    #     lensA = splitsA.apply(pd.Series).stack().str.len()
    #     lensB = splitsB.apply(pd.Series).stack().str.len()
    #     maxA = lensA.max()
    #     maxB = lensB.max()
    #     self.consecutive_nucleations = max(maxA, maxB)
    #     return self.consecutive_nucleations, self.total_nucleations   
