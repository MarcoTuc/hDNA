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
    def __init__(self, structure):
        self.str = structure
        self.left, self.right = structure.split('+')
        self.length = len(self.left)
        if self.length != len(self.right):
            raise BrokenPipeError('Left and Right strands should have the same length')
        self.totbp = sum([True if i != '.' else False for i in self.left])
        if self.totbp != sum([True if i != '.' else False for i in self.right]):
            raise BrokenPipeError('Left and Right base pairs should always match')
        
        self.lì = ''.join(['ì' if i == '(' else '.' for i in self.left])
        self.rì = ''.join(['ì' if i == ')' else '.' for i in self.right])
        self.ì = '+'.join([self.lì, self.rì])
        self.register = 0 if self.duplex else self.get_register()
        self.tail = self.length - abs(self.register)
    
    def get_register(self):
        def shift(s, d):
            if d == +1:
                if s[-1] == 'ì':
                    return s
                else:
                    shf = s[len(s)-1:] + s[:len(s)-1]
                    return shf
            if d == -1:
                if s[0] == 'ì':
                    return s
                else:
                    shf = s[1:] + s[:1] 
                    return shf
        il = 0; ir = 0; dl = None; dr = None
        bb1 = self.rì[::-1]
        bb2 = self.rì[::-1]
        while bb1 != self.lì:
            il += 1
            nbb1 = shift(bb1, 1)
            if nbb1 == self.lì: 
                dl = -il
                break
            if bb1 == nbb1: 
                break
            else: bb1 = nbb1
        while bb2 != self.lì:
            ir += 1
            nbb2 = shift(bb2, -1)
            if nbb2 == self.lì:
                dr = ir 
                break
            if bb2 == nbb2: 
                break
            else: bb2 = nbb2
        dist = dl if dl != None else dr 
        return dist
    
    @property
    def duplex(self):
        L = [True if i != '.' else False for i in self.left]
        R = [True if i != '.' else False for i in self.right]
        if all(L) and all(R):
            return True 