import re

from .strand import Strand
from .complex import Complex
from .model import Model

class Chamber(object):

    """ The chamber object represents all the possible
    nucleation states for the two given single strands """

    def __init__(self, model: Model, s1: Strand, s2: Strand):
        
        self.model = model
        self.s1 = s1
        self.s2 = s2.invert
        self.duplex = Complex(self.model, self.s1, self.s2, duplex=True)
        """ General Slidings """
        self.compute_slidings_structured(3)
        """ Off-Register Nucleation Cores """
        self.compute_offcores(3)


    def compute_offcores(self, min_nucleation):
        self.offcores = []
        for complex in self.slidings:
            if complex.consecutive_nucleations >= min_nucleation:
                self.offcores.append(complex)
    
    def compute_oncores(self, min_nucleation):
        self.oncores = []
        self.compute_slidings_structured(min_nucleation)
        for complex in self.slidings:
            if complex.consecutive_nucleations >= min_nucleation:
                self.oncores.append(complex)
    
    def compute_slidings_structured(self, min_nucleation):
        n = min_nucleation
        self.slidings = []
        for b in range(n,self.s1.length): 
            slidingstruct = "i"*b+"."*(self.s1.length- b)
            slidingstruct = slidingstruct+"+"+slidingstruct
            structureout = self.parse_structure(slidingstruct, self.s1, self.s2)
            self.slidings.append(Complex(self.model, self.s1, self.s2, structure=structureout, offregister=True))
        fullstructure = "("*self.s1.length+"+"+")"*self.s2.length
        self.slidings.append(Complex(self.model, self.s1, self.s2, structure=fullstructure, duplex=True))
        for b in range(1, self.s1.length - n + 1):
            slidingstruct = "."*b+"i"*(self.s1.length-b)
            slidingstruct = slidingstruct+"+"+slidingstruct
            structureout = self.parse_structure(slidingstruct, self.s1, self.s2)
            self.slidings.append(Complex(self.model, self.s1, self.s2, structure=structureout, offregister=True))
    
    
    def parse_structure(self, structure, seq_a, seq_b):
        struct_a, struct_b = structure.split('+')
        seq_a = seq_a.sequence
        seq_b = seq_b.sequence 
        sx, dx = self.structurecut(seq_a, seq_b[::-1], struct_a, struct_b[::-1])
        patch_sx = ''.join(['i' for i in range(len(sx))])
        patch_dx = ''.join(['i' for i in range(len(dx))])
        out1 = re.sub(patch_sx,sx,struct_a)
        out2 = re.sub(patch_dx,dx,struct_b)
        structureout = out1+"+"+out2
        return structureout

    def structurecut(self, string1, string2, structure1, structure2):
        cut1 = ''.join([n1 for n1, s1 in zip(string1, structure1) if s1 == 'i'])
        cut2 = ''.join([n2 for n2, s2 in zip(string2, structure2) if s2 == 'i'])
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


    def iswattsoncrick(self,n1,n2):
        if (n1 == 'A' and n2 == 'T') or (n1 == 'C' and n2 == 'G') or (n1 == 'T' and n2 == 'A') or (n1 == 'G' and n2 == 'C'):
            return True
        else:
            return False 
            
    def convolve(self, min_nucleation):
        """Return the onregister nucleations"""
        pass

    def duplexenergy(self):
        return self.duplex.fenergy