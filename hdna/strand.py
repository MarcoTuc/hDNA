from random import choice 
from .model import Model

class Strand(object):

    def __init__(self, model: Model, sequence: str, direction='53', bottom=None):

        #TODO: Check if the sequence has non WT stuff 

        if type(model) != Model:
            raise TypeError("Model must be an instance of hdna.Model")
        self.model = model
        
        if type(sequence) != str:
            raise TypeError("Sequence must be a string")
        self.sequence = sequence
        self.length = len(sequence)
        self.direction = direction
        self.sdist = [i for i in range(self.length)]
        if bottom == self.direction[1]:
            self.sdist = self.sdist[::-1]

    def dimension(self):
        """ I don't remember what I wanted to do here"""
        pass

    def nucleotides(self):
        pass
    
    def cut(self, start, stop):
        return Strand(self.model, self.sequence[start:stop])

    def complementary(self, **kwargs):
        wc = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        seq = ''.join([wc[self.sequence[i]] for i in range(self.length)])
        dc = {None: None, '53':'35', '35':'53'}
        return Strand(self.model, seq, direction=dc[self.direction], **kwargs)
    
    def random(model, length, direction='53', **kwargs):
        """ Generate a random sequence of given length """
        seq = ''.join([choice(['A','T','C','G']) for i in range(length)])
        return Strand(model, seq, direction, **kwargs)

    def get_ix(string, char):
        indices = []
        for i, e in enumerate(list(string)):
            if e == char:
                indices.append(i)
        return indices 
    
    @property
    def invert(self):
        return Strand(self.model, self.sequence[::-1])
    


class Structure(object):
    """ If I have time I will put here all the structure related Methods
        that for now are just copy-pasted around all classes """
    pass 