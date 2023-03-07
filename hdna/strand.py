import numpy as np 

from random import choice 
from .model import Model
from .params import * 

class Strand(object):

    def __init__(self, model: Model, sequence: str, direction=None):

        #TODO: Check if the sequence has non WT stuff 
        if type(model) != Model:
            raise TypeError("Model must be an instance of hdna.Model")
        self.model = model
        
        if type(sequence) != str:
            raise TypeError("Sequence must be a string")
        self.sequence = sequence
        
        if direction is not None: 
            if self.direction == '53':
                self.fivethree = self.sequence
                self.threefive = self.invert.sequence
            
            if self.direction == '35':
                self.fivethree = self.invert.sequence
                self.threefive = self.sequence 

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

    def complementary(self, direction=None):
        wc = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        seq = ''.join([wc[self.sequence[i]] for i in range(self.length)])
        dc = {None: None, '53':'35', '35':'53'}
        return Strand(self.model, seq, dc[direction])
    
    def random(model, length, direction=None):
        """ Generate a random sequence of given length """
        seq = ''.join([choice(['A','T','C','G']) for i in range(length)])
        return Strand(model, seq, direction)

    @property
    def invert(self):
        return Strand(self.model, self.sequence[::-1])

   


class Structure(object):
    def __init__(self, structure, fromtable=False):
        if fromtable:
            self.table = structure
            self.left = ''.join(['(' if e else '.' for e in structure[:len(structure)//2]])
            self.right = ''.join([')' if e else '.' for e in structure[1+len(structure)//2:]])
            self.str = '+'.join([self.left, self.right])
        else:
            self.str = structure
            self.left, self.right = structure.split('+')
            self.table  = [True if i not in ['.','+'] else False for i in self.str]
        
        self.length = len(self.left)
        if self.length != len(self.right):
            raise BrokenPipeError(f'Left and Right strands should have the same length: {self.str}')
        self.totbp = sum(self.table)/2
        if self.totbp != sum([True if i != '.' else False for i in self.right]):
            raise BrokenPipeError('Left and Right base pairs should always match')
        
        if any(self.table):
            self.lì = ''.join(['ì' if i == '(' else '.' for i in self.left])
            self.rì = ''.join(['ì' if i == ')' else '.' for i in self.right])
            self.ì = '+'.join([self.lì, self.rì])
            
            self.register = 0 if self.duplex else self.get_register()        
            self.get_geometry()
            self.get_pktails()
            self.maxtails()
            self.overlappingspheres()
            self.inchwormingtails()
        else:
            self.ss = True

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
    
    def get_geometry(self):
        tail_ll = 0
        tail_lr = 0
        i = 0
        ell = self.left[i]
        while ell == '.':
            i += 1
            ell = self.left[i]
            tail_ll += 1
        i = 1
        elr = self.left[-i]
        while elr == '.':
            i += 1
            elr = self.left[-i]
            tail_lr += 1
        bulkl = self.length - tail_ll - tail_lr

        tail_rl = 0
        tail_rr = 0
        i = 0
        erl = self.right[i]
        while erl == '.':
            i += 1
            erl = self.right[i]
            tail_rl += 1
        i = 1
        err = self.right[-i]
        while err == '.':
            i += 1
            err = self.right[-i]
            tail_rr += 1
        bulkr = self.length - tail_rl - tail_rr

        if bulkl != bulkr:    
            self.geometry = {'left': (tail_ll, bulkl, tail_lr),
                            'right':(tail_rl, bulkr, tail_rr)}
        else: 
            
            self.tail_ll = tail_ll
            self.tail_lr = tail_lr
            self.tail_rl = tail_rl
            self.tail_rr = tail_rr
            
            self.tails   = {'ll': self.tail_ll, 'lr': self.tail_lr,
                            'rl': self.tail_rl, 'rr': self.tail_rr}
            
            self.tails_l = {'l': self.tail_ll, 'r': self.tail_lr}
            self.tails_r = {'l': self.tail_rl, 'r': self.tail_rr}

            self.bulk    = bulkl
            self.geometry = {'register': self.register,
                             'left': self.tails_l,
                             'bulk': self.bulk,
                             'right': self.tails_r}

    def get_pktails(self):
        if self.register > 0:
            tls = 'r'
            self.pktail_l = self.tails_l[tls]
            self.pktail_r = self.tails_r[tls]
        elif self.register < 0:
            tls = 'l'
            self.pktail_l = self.tails_l[tls]
            self.pktail_r = self.tails_r[tls]
        else:
            self.pktail_l = 0
            self.pktail_r = 0

    def maxtails(self):
        # for tail in [self.tail_ll, self.tail_lr, self.tail_rl, self.tail_rr]:
        if self.tail_ll > self.tail_lr:
            self.maxtail_l = 'll'
        elif self.tail_ll < self.tail_lr:
            self.maxtail_l = 'lr'
        else: 
            self.maxtail_l = 'both'
        
        if self.tail_rl > self.tail_rr:
            self.maxtail_r = 'rl'
        elif self.tail_rl < self.tail_rr:
            self.maxtail_r = 'rr'
        else: 
            self.maxtail_r = 'both'
        
    def sumtails(self):
        return sum([self.tail_ll, self.tail_lr, self.tail_rl, self.tail_rr])

    def spheresoverlap(self, r1, r2, d):
        separation = r1+r2-d
        if separation <= 0:
            return 0 
        else: 
            term1 = separation**2
            term2 = (d**2 + 2*d*r1 - 3*(r1**2) + 2*d*r2 + 6*r1*r2 - 3*(r2**2))
            overlap = np.pi*term1*term2/(12*d)
            sphere1  = (4/3)*np.pi*(r1**3)
            sphere2  = (4/3)*np.pi*(r2**3)
            normoverlap = overlap/(sphere1+sphere2)
            return normoverlap 

    def overlappingspheres(self):
        d  = self.bulk * DXGEO.MONODIST #stiff rod bro (persistence duplex >> persistence simplex)
        if self.tail_ll > 0 and self.tail_rl > 0:
            # left up and right low spheres 
            r1 = np.sqrt(self.tail_ll*(SXGEO.MONODIST**2))
            r2 = np.sqrt(self.tail_rl*(SXGEO.MONODIST**2))
            overlap_lurd = self.spheresoverlap(r1, r2, d)
        else:
            overlap_lurd = 0
        if self.tail_rr > 0 and self.tail_lr > 0:
            # left low and right up spheres 
            r3 = np.sqrt(self.tail_rr*(SXGEO.MONODIST**2))
            r4 = np.sqrt(self.tail_lr*(SXGEO.MONODIST**2))
            overlap_ldru = self.spheresoverlap(r3, r4, d)    
        else:
            overlap_ldru = 0
        overlap = (overlap_lurd + overlap_ldru)/2
        self.pkoverlap = overlap 

    def inchwormingtails(self):
        if self.tail_ll != 0 and self.tail_rr != 0:
            iw_left = (self.tail_ll+self.tail_rr)/(2*self.length*abs(self.register))  #0.5*abs(self.tail_ll+self.tail_rr)/abs(self.register)
        else: iw_left = 0
        if self.tail_lr != 0 and self.tail_rl != 0:
            iw_right = (self.tail_lr+self.tail_rl)/(2*self.length*abs(self.register)) #0.5*abs(self.tail_lr+self.tail_rl)/abs(self.register)
        else: iw_right = 0 
        self.iw_left = iw_left          #abs(self.tail_ll+self.tail_rr)/iw_left #what about a product here
        self.iw_right = iw_right        #abs(self.tail_lr+self.tail_rl)/iw_right #what about a product here
        return self.iw_left, self.iw_right
    
    @property
    def inchwormingbulge(self):
        return abs(1/self.register)

    @staticmethod
    def empty(length):
        return Structure('.'*length+'+'+'.'*length)  

    @property
    def duplex(self):
        L = [True if i != '.' else False for i in self.left]
        R = [True if i != '.' else False for i in self.right]
        if all(L) and all(R):
            return True 
        