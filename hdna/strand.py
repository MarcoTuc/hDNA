from numpy.random import choice 
from .model import Model

class Strand(object):

    def __init__(self, model: Model, sequence: str):

        #TODO: Check if the sequence has non WT stuff 
        if type(model) != Model:
            raise TypeError("Model must be an instance of hdna.Model")
        self.model = model
        
        if type(sequence) != str:
            raise TypeError("Sequence must be a string")
        self.sequence = sequence

        self.length = len(sequence)

    def secstruct(self):
        """ method for getting the most expressed 
        secondary structures. These structures will be objects
        themselves with their own shit """       
        pass

    def dimension(self):
        """ I don't remember what I wanted to do here"""
        pass

    def nucleotides(self):
        pass
    
    def cut(self, start, stop):
        return Strand(self.model, self.sequence[start:stop])

    def complementary(self):
        wc = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        seq = ''.join([wc[self.sequence[i]] for i in range(self.length)])
        return Strand(self.model, seq)
    
    def random(model, length):
        """ Generate a random sequence of given length """
        seq = ''.join([choice(['A','T','C','G']) for i in range(length)])
        return Strand(model, seq)


    @property
    def invert(self):
        return Strand(self.model, self.sequence[::-1])


class Structure(object):
    """ If I have time I will put here all the structure related Methods
        that for now are just copy-pasted around all classes """
    pass 