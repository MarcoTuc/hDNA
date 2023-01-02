import nupack as nu 
import classes.strand as ss
import classes.model as m

class Complex(object):
    
    def __init__(   self, 
                    model: m.Model, 
                    s1: ss.Strand, 
                    s2: ss.Strand, 
                    structure=None,
                    duplex=False,
                    offregister=False,
                    onregister=False,
                    dangsx=False,
                    dangdx=False
                ):
        
        self.model = model

        self.s1 = s1
        self.s2 = s2 

        self.l1 = s1.length
        self.l2 = s2.length

        self.mismatches = []
        self._get_mismatches()

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
        self.fenergy = j[1]
        self.dpfunc = float(j[0])

        """ These variables are used to say if the current
        Coils object is duplexed, is offregister or is onregister"""

        self.danglesx = dangdx
        self.dangledx = dangsx

        self.structure = structure

        self.maxnucleation()

        self.duplex = duplex
        self.offregister = offregister
        self.onregister = onregister


    def _get_mismatches(self):
        for nuc1, nuc2 in zip(self.s1.sequence, self.s2.sequence):
            self.mismatches.append(self._iswattsoncrick(nuc1, nuc2))

    def _iswattsoncrick(self,n1,n2):
        """
        check for mismatches
        """
        if (n1 == 'A' and n2 == 'T') or (n1 == 'C' and n2 == 'G') or (n1 == 'T' and n2 == 'A') or (n1 == 'G' and n2 == 'C'):
            return True
        else:
            return False 
    
    def gotnucleationhere(self, mincore):
        self.nucleation = None
        if self.l1 != self.l2:
            return ValueError('This method cannot be used for complexes made with different length strands')
        if self.nucleation_size >= mincore:
            self.nucleation == True
        else:
            self.nucleation = False

    def maxnucleation(self):
        self.nucleation_size = 0
        for n in range(0, len(self.mismatches) + 1):
            if self.checknucleation(n):
                self.nucleation_size = n
        return self.nucleation_size

    def checknucleation(self, mincore):
        self.nucleation = None
        if self.l1 != self.l2:
            return ValueError('This method cannot be used for complexes made with different length strands')
        for b in range(len(self.mismatches) - (mincore-1)):
            if all(self.mismatches[b:b+mincore]):
                return True 
        return False 
    
    def correctedstrand(self):
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

    @property
    def sequences(self):
        return self.s1.sequence+"+"+self.s2.sequence

