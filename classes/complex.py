import nupack as nu 
import numpy as np 
import pandas as pd 

from classes.strand import Strand
from classes.model import Model

class Complex(object):
    
    def __init__(   self, 
                    model: Model, 
                    s1: Strand, 
                    s2: Strand, 
                    structure=None,
                    duplex=False,
                    offregister=False,
                    onregister=False,
                ):
        
        self.model = model

        self.s1 = s1
        self.s2 = s2 

        self.l1 = s1.length
        self.l2 = s2.length

        self.mismatches = []
        self._get_mismatches()

        """ These variables are used to say if the current
        Coils object is duplexed, is offregister or is onregister"""

        self.duplex = duplex
        self.offregister = offregister
        self.onregister = onregister

        self.structure = structure

        self.getnupackproperties()

        if self.duplex == True:
            #TODO update this when considering mismatches 
            self.consecutive_nucleations = min(self.l1, self.l2)
            self.total_nucleations = self.consecutive_nucleations
        else:
            self.nucleationsize()
    
    def getnupackproperties(self):
            """ Nupack related properties """
            self.nuStrand1 = nu.Strand(self.s1.sequence, name = 'a')
            self.nuStrand2 = nu.Strand(self.s2.invert.sequence, name = 'b') # An inversion here is needed because in this program strands are defined as 5-3 against 3-5 but in NUPACK all strands are defined 5-3 and the program takes care to turn them around and so on
            self.nuComplex = nu.Complex([self.nuStrand1,self.nuStrand2], name = 'c')

            """For now I will keep NUPACK ensemble model variable
            fixed at the following value:
            ensemble = 'stacking' """

            nupackmodel = nu.Model(material=self.model.material, 
                            ensemble='stacking', 
                            celsius=self.model.celsius, 
                            sodium=self.model.Na, 
                            magnesium= self.model.Mg)

            j = nu.pfunc(self.nuComplex, nupackmodel)
            self.G = j[1]
            self.Z = float(j[0])

    def nucleationsize(self):
        """ This method gives back two quantities:
            - Total number of non-mismatched nucleated base pairs
            - Size of the biggest consecutive sequence of nucleated base pairs """

        A, B = self.splitstructure()
        A, B = pd.Series(A), pd.Series(B)

        nA = A.str.count("\(").sum()
        nB = B.str.count("\)").sum()

        if nA != nB: raise ValueError('total number of based pairs should match for each strand')
        else: self.total_nucleations = nA 

        splitsA = A.str.split('[^(]')
        splitsB = B.str.split('[^)]')

        lensA = splitsA.apply(pd.Series).stack().str.len()
        lensB = splitsB.apply(pd.Series).stack().str.len()

        maxA = lensA.max()
        maxB = lensB.max()

        self.consecutive_nucleations = max(maxA, maxB)

        return self.consecutive_nucleations, self.total_nucleations
    

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

    def splitstructure(self):
        self.struct_a, self.struct_b = self.structure.split('+')
        return self.struct_a, self.struct_b

    @property
    def sequences(self):
        return self.s1.sequence+"+"+self.s2.sequence


  
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