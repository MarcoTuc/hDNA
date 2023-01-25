import re

from .strand import Strand
from .complex import Complex, Sliding 
from .model import Model


class Chamber(object):

    """ The chamber object contains every possible
    nucleation state for the two given single strands """

    def __init__(self, model: Model, s1: Strand, s2: Strand, mincore):

        self.model = model
        self.s1 = s1
        self.s2 = s2.invert

        self.singlestranded = Complex(self.model, self.s1, self.s2, singlestranded=True)
        self.mincore = mincore
        #TODO: generalize the minimum nucleation size, right now it is just 3

        """ General Slidings """
        self.compute_slidings_structured()
        
        """ Off-Register Nucleation Cores """
        self.compute_offcores()

        """ On-Register Nucleation Cores """
        self.compute_oncores()

        """
        self.finalstructure --->    compute the most stable structure for the 
                                    given sequence, this will be the arrival 
                                    of trajectories on the kinetic simulation
        """



#########################################
##### Non-Native Nucleation Methods #####
#########################################

    def compute_offcores(self):
        self.offcores = []
        for complex in self.slidings:
            if complex.consecutive_nucleations >= self.mincore:
                self.offcores.append(complex)

    def compute_slidings_structured(self, verbose=False):
        n = self.mincore
        self.slidings = []
        for b in range(n, self.s1.length): 
            slidingstruct = "ì"*b+"."*(self.s1.length-b)
            slidingstruct = slidingstruct + "+" + slidingstruct
            structureout = self.parse_structure(slidingstruct, self.s1, self.s2)
            if verbose: print(self.s1.sequence+'+'+self.s2.sequence); print(slidingstruct, (self.s1.length - b)); print(structureout,'\n')
            self.slidings.append(Sliding(self.model, self.s1, self.s2, structure=structureout, offregister="left", dpxdist=(self.s1.length - b)))
        fullstructure = "("*self.s1.length+"+"+")"*self.s2.length
        if verbose: print(fullstructure)
        self.duplex = Complex(self.model, self.s1, self.s2, structure=fullstructure, duplex=True)
        for b in range(1, self.s1.length - n + 1):
            slidingstruct = "."*b+"ì"*(self.s1.length-b)
            slidingstruct = slidingstruct+"+"+slidingstruct
            structureout = self.parse_structure(slidingstruct, self.s1, self.s2)
            if verbose: print(self.s1.sequence+'+'+self.s2.sequence); print(slidingstruct, b); print(structureout,'\n')
            self.slidings.append(Sliding(self.model, self.s1, self.s2, structure=structureout, offregister='right', dpxdist=b))

    def split_offcores(self):
        n = int(len(self.offcores)/2)
        return self.offcores[:n], self.offcores[n:][::-1]


#####################################
##### Native Nucleation Methods #####
#####################################

    def compute_oncores(self):
        self.oncores = []
        self.native_nucleation_structures()
        for structureì in self.nativeì:
            structureout = self.parse_structure(structureì, self.s1, self.s2)
            self.oncores.append(Complex(self.model, self.s1, self.s2, structure=structureout, onregister=True))
        self.oncores = [core for core in self.oncores if core.total_nucleations >= self.mincore]
        return self.oncores

    def clean_oncores(self):
        for core in self.oncores:
            if core.total_nucleations < self.mincore:
                self.oncores.pop(core)

    def native_nucleation_structures(self):
        """Return the onregister nucleations"""
        n = self.mincore
        self.nativeì = []
        for i in range(self.s1.length - n + 1):
            nucleation = "".join(["." * i, "ì" * n,"." * (self.s1.length - i - n), "+", "." * (self.s2.length - i - n), "ì" * n, "." * i])
            self.nativeì.append(nucleation)
        return self.nativeì


#############################################
##### General Structure Parsing Methods #####
#############################################
    
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
        patch_sx = ''.join(['ì' for i in range(len(sx))])
        patch_dx = ''.join(['ì' for i in range(len(dx))])
        out1 = re.sub(patch_sx,sx,struct_a)
        out2 = re.sub(patch_dx,dx,struct_b)
        structureout = out1+"+"+out2
        return structureout


    def structurecut(self, string1, string2, structure1, structure2):
        cut1 = ''.join([n1 for n1, s1 in zip(string1, structure1) if s1 == 'ì'])
        cut2 = ''.join([n2 for n2, s2 in zip(string2, structure2) if s2 == 'ì'])
        sx = ''
        dx = ''
        # print(cut1, cut2)
        for n1, n2 in zip(cut1, cut2):
            if self.iswattsoncrick(n1, n2):
                sx += '('
                dx += ')'
            else: 
                sx += '.'
                dx += '.' 
        return sx, dx[::-1]


##########################
##### Helper Methods #####
##########################

    def iswattsoncrick(self,n1,n2):
        if (n1 == 'A' and n2 == 'T') or (n1 == 'C' and n2 == 'G') or (n1 == 'T' and n2 == 'A') or (n1 == 'G' and n2 == 'C'):
            return True
        else:
            return False 

    def duplexenergy(self):
        return self.duplex.fenergy