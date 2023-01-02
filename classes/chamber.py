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
        self.compute_slidings(3)

        """ Off-Register Nucleation Cores """
        self.compute_offcores(3)


    def compute_offcores(self, min_nucleation):
        n = min_nucleation
        self.offcores = []
        self.compute_slidings(min_nucleation)
        for complex in self.slidings:
            if complex.nucleation_size >= min_nucleation:
                self.offcores.append(complex)
    
    def compute_oncores(self, min_nucleation):
        n = min_nucleation
        self.oncores = []
        self.compute_slidings(min_nucleation)
        for complex in self.slidings:
            if complex.nucleation_size >= min_nucleation:
                self.oncores.append(complex)


    def compute_slidings(self, min_nucleation):
        n = min_nucleation
        self.slidings = []
        for b in range(n,self.s1.length): 
            seq1 = self.s1.cut(0,b)
            seq2 = self.s2.invert.cut(0,b).invert
            self.slidings.append(Complex(self.model, seq1, seq2, offregister=True))
        self.slidings.append(Complex(self.model, self.s1, self.s2, duplex=True))
        for b in range(1, self.s1.length - n + 1):
            seq1 = self.s1.cut(b,None)
            seq2 = self.s2.cut(0,self.s2.length - b)
            self.slidings.append(Complex(self.model, seq1, seq2, offregister=True))
    
    def compute_slidings_structured(self, min_nucleation):
        n = min_nucleation
        self.slidings = []
        for b in range(n,self.s1.length): 
            structure = "^"*b+"."*(self.s1.length- b)
            structure = structure+"+"+structure
            # structure = self.parse_structure(structure, self.s1, self.s2)
            self.slidings.append(Complex(self.model, self.s1, self.s2, structure, offregister=True))
        fullstructure = "("*self.s1.length+"+"+")"*self.s2.length
        self.slidings.append(Complex(self.model, self.s1, self.s2, fullstructure, duplex=True))
        for b in range(1, self.s1.length - n + 1):
            structure = "."*b+"^"*(self.s1.length-b)
            structure = structure+"+"+structure
            # structure = self.parse_structure(structure, self.s1, self.s2)
            self.slidings.append(Complex(self.model, self.s1, self.s2, structure, offregister=True))
    
    
    def parse_structure(self, structure, seq_a, seq_b):

        struct_a, struct_b = structure.split('+')

        seq_a = seq_a.sequence
        seq_b = seq_b.sequence 

        dotbracket_a = "" #left dotbracket
        dotbracket_b = "" #right dotbracket 

        overlap1, overlap2 = self.structurecut(seq_a, seq_b, struct_a, struct_b)

        pairs = zip(seq_a, seq_b, struct_a, struct_b)
        for i, (base_a, base_b, struct_a, struct_b) in enumerate(pairs):
            if struct_a == ".": dotbracket_a += '.'
            if struct_b == ".": dotbracket_b += '.'
            if struct_a == "^":
                if self.iswattsoncrick(base_a,base_b): dotbracket_a += "("
                else: dotbracket_a += "."
            if struct_b == "^":
                if self.iswattsoncrick(base_a,base_b): dotbracket_b += ")"
                else: dotbracket_b += "."
        return dotbracket_a+'+'+dotbracket_b

    def structurecut(self, string1, string2, structure1, structure2):
        cut1 = ''.join([s1 for s1, s2 in zip(string1, structure1) if s2 == '^'])
        cut2 = ''.join([s2 for s1, s2 in zip(structure2, string2) if s1 == '^'])
        return cut1, cut2[::-1]
                
    def iswattsoncrick(self,n1,n2):
        """
        check for mismatches
        """
        if (n1 == 'A' and n2 == 'T') or (n1 == 'C' and n2 == 'G') or (n1 == 'T' and n2 == 'A') or (n1 == 'G' and n2 == 'C'):
            return True
        else:
            return False 
            

    def convolve(self, min_nucleation):
        """Return the onregister nucleations"""
        pass

    def duplexenergy(self):
        return self.duplex.fenergy