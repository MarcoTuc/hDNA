class Strand(object):

    def __init__(self, model, sequence):

        #TODO: Check if the sequence has non WT stuff 
        self.model = model
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

    @property
    def invert(self):
        return Strand(self.model, self.sequence[::-1])

