
import sys
sys.path.insert(0, '.')

from values import *
from kineticsfunctions import *
from conf import *

TEMPERATURE = 273.15+26


def iswc(n1,n2):
    """
    check for mismatches
    """
    if n1 == 'A' and n2 == 'T' or n1 == 'C' and n2 == 'G' or n1 == 'T' and n2 == 'A' or n1 == 'G' and n2 == 'C':
        return True
    else:
        return False


def posab(n1,n2):
    """
    check for abasic damages
    """
    if n1 == 'U':
        return int(0)
    elif n2 == 'U':
        return int(1)
    else:
        return int

def wc(base:str):
    """
    return the corresponding Watson-Crick opposite base
    """
    if base == 'A': return str('T')
    if base == 'T': return str('A')
    if base == 'C': return str('G')
    if base == 'G': return str('C')

def complementary_strand(strand:str):
    return ''.join([wc(strand[i]) for i in range(len(strand))])

def modified_matrix(strandA, strandB, verbose=True):
    """
    when calculating free energies we first compute them
    as the strands are perfectly watson-crick and then correct their 
    free energy by accounting for damages. 
    This function returns a corrected WC strand for performing
    free energy calculations using the nearest-neighbor model from SantaLucia
    """
    strands = [list(strandA),list(strandB)]
    if verbose == True:
        print('The original strand matrix is \n',strands[0],'\n',strands[1])
    for i in [0,1]:
        for n in range(len(strands[i])):
            if strands[i][n] == 'U':
                strands[i][n] = wc(strands[i-1][n])
    if verbose == True:
        print('The corrected strand matrix is \n',strands[0],'\n',strands[1])
    return strands

def get_flag(strandA, strandB):
    """
    flags are used to tell the algorithm if the strand presents any 
    mismatches, abasic defects, terminal AT pairs or final abasic sites
    """
    strands = [list(strandA),list(strandB)]
    flag = {'abasic':[],'abasic_final':False,'mismatches':[],'termAT':False}
    for i in [0,1]:
        for n in range(1,len(strands[i])-1):
            if strands[i][n] == 'U':
                flag['abasic'].append([i,n])
        for n in range(len(strands[i])):
            if strands[i][n] != wc(strands[i-1][n]):
                flag['mismatches'].append([i,n])
            if strands[i][0] or strands[i][-1] == 'A' or 'T':
                flag['termAT'] = True
    return flag


def abasic_ddg(context,opposing_nucleotide):
    """
    compute the free energy penalty for abasic sites based on the 
    sequence context surrounding the abasic site
    """
    if context in [['G','U','G'],['T','U','G'],['G','U','T'],['T','U','T']]:
        return ab['GUG'][opposing_nucleotide]
    if context in [['C','U','C'],['A','U','C'],['C','U','A'],['A','U','A']]:
        return ab['CUC'][opposing_nucleotide]
    else:
        avg = (ab['GUG'][opposing_nucleotide]+ab['CUC'][opposing_nucleotide])/2
        return avg

def add_initiation(free_energy_list,temperature,verbose=True):
    """
    just adding the initiation energy needed for computing the total free energy
    """
    free_energy_list.append(gibbs_free(hs['initiation']['H'],hs['initiation']['S'],temperature=temperature))
    if verbose == True:
        print('Adding initiation free energy equal to', gibbs_free(hs['initiation']['H'],hs['initiation']['S'],temperature=temperature))


def terminalATpenalty(strands, free_energy_list, temperature,verbose=True):
    """
    adding the terminal AT penalty if the terminal AT penalty flag is true
    """
    terminal = gibbs_free(hs['TerminalATPenalty']['H'],hs['TerminalATPenalty']['S'],temperature=temperature)
    if strands[0][0] == 'A' and strands[1][0] == 'T':
        free_energy_list.append(terminal)
        if verbose == True:
            print('Adding terminal AT penalty for initial bases equal to:',terminal)
    if strands[0][0] == 'T' and strands[1][0] == 'A':
        free_energy_list.append(terminal)
        if verbose == True:
            print('Adding terminal AT penalty for initial bases equal to:',terminal)
    if strands[0][-1] == 'A' and strands[1][-1] == 'T':
        free_energy_list.append(terminal)
        if verbose == True:
            print('Adding terminal AT penalty for final bases equal to:',terminal)
    if strands[0][-1] == 'T' and strands[1][-1] == 'A':
        free_energy_list.append(terminal)
        if verbose == True:
            print('Adding terminal AT penalty for final bases equal to:',terminal)


##### MAIN FUNCTION 

def string_free(strandA, strandB, temperature=TEMPERATURE,verbose=False):
    """
    This is the main function for calculating the whole strand free energy
    given the initial input of the two strands we want to hybridize. 
    """
    if verbose == True:
        print('All values are measured in Kcal/mol')
    #print('Dangling end contribution is approximated as an additional AT basepair at the end: https://academic.oup.com/nar/article/28/9/1929/2903888')
    strands = modified_matrix(strandA,strandB,verbose=verbose)
    original_strands = [list(strandA),list(strandB)]
    flag = get_flag(strandA, strandB)
    if verbose == True:
        print(flag)
    pairs_dg = []
    abasic_penalties = []
    add_initiation(pairs_dg,temperature,verbose=verbose)
    if flag['termAT'] == True:
        terminalATpenalty(strands, pairs_dg, temperature,verbose=verbose)
    for i in range(int(len(strands[0]))-1):
        if strands[0][i] == wc(strands[1][i]):
            pairfree = gibbs_free(hs[strands[0][i]][strands[0][i+1]][strands[1][i+1]]['H'],hs[strands[0][i]][strands[0][i+1]][strands[1][i+1]]['S'],temperature)
            pairs_dg.append(pairfree)
            if verbose == True:
                if strands[0][i+1] == wc(strands[1][i+1]):
                    print('Delta G for the WC pair',strands[0][i],strands[0][i+1],'-',strands[1][i],strands[1][i+1],':',pairfree)
                else:
                    print('Delta G for the mismatch',strands[0][i],strands[0][i+1],'-',strands[1][i],strands[1][i+1],':',pairfree)
        else:
            mismfree = gibbs_free(hs[strands[1][i+1]][strands[1][i]][strands[0][i]]['H'],hs[strands[1][i+1]][strands[1][i]][strands[0][i]]['S'],temperature)
            pairs_dg.append(mismfree)
            if verbose == True:
                print('Delta G for the mismatch',strands[0][i],strands[0][i+1],'-',strands[1][i],strands[1][i+1],':',mismfree)
    if flag['abasic'] != []:
        for ab in flag['abasic']:
            context = [original_strands[ab[0]][ab[1]+j] for j in [-1,0,1]]
            opp_base = original_strands[ab[0]-1][ab[1]]
            abasic_penalties.append(-abasic_ddg(context,opp_base))
            if verbose == True:
                print('Abasic defect found in strand',ab[0],'at position',ab[1])
                print('Context:',''.join(context)+'; Opposing base:',opp_base)
                print('The free energy penalty is',-abasic_ddg(context,opp_base))
    if 'U' in [original_strands[0][0],original_strands[1][0]]:
        abasic_penalties.append(2)
        if verbose == True:
            print('Added a dangling abasic penalty of 2 at initial position',original_strands[0][0] or original_strands[1][0])
    if 'U' in [original_strands[0][-1],original_strands[1][-1]]:
        abasic_penalties.append(2)
        if verbose == True:
            print('Added a dangling abasic penalty of 2 at final position')

    repaires_freene = sum(pairs_dg)
    tot_abasic_penalty = sum(abasic_penalties)
    freene = repaires_freene + tot_abasic_penalty
    if verbose == True:
        print('Delta G of repaired strand:',repaires_freene)
        print('Total abasic penalty:',tot_abasic_penalty)
        print('Final Delta G:',freene)
    return freene

class BreakPointError(Exception):
    def __init__(self, message='Break point cannot be larger than string'):
        self.message = message
        super().__init__(self.message)
    def __str__(self) -> str:
        return f'{self.message}'
    pass


############## LIDA Kinetics function

def get_lida_rates_fixedforward(Tl,Or1,Or2,breaking_point,invert=False,temperature=TEMPERATURE,verbose=True):
    """
    Here all the kinetic rates are computed as a function of the template
    and oligomer strands given for a LIDA system
    """
    if breaking_point > len(Tl):
        raise BreakPointError

    kf = 2e7
    
    Tr = list(''.join([Or1,Or2]))
    Tl = list(Tl)                       # 5'-3'   NNNNNNNNNNNNNNNNNNNNNN
    Or1 = list(Or1)                     # 3'-5'   NNNNNNNNNN
    Or2 = list(Or2)                     # 3'-5'             NNNNNNNNNNNN
    if invert == True:
        Or1 = list(Or1)[::-1]       
        Or2 = list(Or2)[::-1]           
    Ol1 = Tl[:breaking_point]                          
    Ol2 = Tl[breaking_point:]       

    space1 = ' '.join(['' for i in range(2*len(Or1))])
    space2 = ' '.join(['' for i in range(2*len(Ol1))])

    if verbose == True:
        print(' D: destabilizing \n'
            ' Template:     ',' '.join(Tl),'\n',
            'D-Oligomer1:  ',' '.join(Or1),'\n',
            'D-Oligomer2:  ',space1,' '.join(Or2),'\n',
            'D-Template:   ',' '.join(Tr),'\n',
            'Oligomer1:    ',' '.join(Ol1),'\n',
            'Oligomer2:    ',space2,' '.join(Ol2),'\n',)

    rates = {   
                'forward' :{},
                'backward':{}
            }
    
    #Routines for getting backward reaction rates

    Tl_dx = Tl[:len(Or1)]
    dg_TlOr1 = string_free(Tl_dx,Or1,temperature,verbose=verbose)
    kd_TlOr1 = k_equilibrium(dg_TlOr1,temperature)
    kb_TlOr1 = kf/kd_TlOr1
    if verbose == True:
        print('TlOr1 Done \n \n')

    Tl_sx = Tl[len(Or1):]
    dg_TlOr2 = string_free(Tl_sx,Or2,temperature,verbose=verbose)
    kd_TlOr2 = k_equilibrium(dg_TlOr2,temperature)
    kb_TlOr2 = kf/kd_TlOr2
    if verbose == True:
        print('TlOr2 Done\n \n')

    Tr_dx = Tr[:len(Ol1)]
    dg_TrOl1 = string_free(Tr_dx,Ol1,temperature,verbose=verbose)
    kd_TrOl1 = k_equilibrium(dg_TrOl1,temperature)
    kb_TrOl1 = kf/kd_TrOl1
    if verbose == True:
        print('TrOl1 Done\n \n')

    Tr_sx = Tr[len(Ol1):]
    dg_TrOl2 = string_free(Tr_sx,Ol2,temperature,verbose=verbose)
    kd_TrOl2 = k_equilibrium(dg_TrOl2,temperature)
    kb_TrOl2 = kf/kd_TrOl2
    if verbose == True:
        print('TrOl2 Done\n \n')

    dg_Td = string_free(Tl,Tr,temperature,verbose=verbose)
    kd_Td = k_equilibrium(dg_Td,temperature)
    kb_Td = kf/kd_Td
    if verbose == True:
        print('Td Done\n \n')
    if verbose == True:
        print('kb_TlOr1:',"{:.3e}".format(kb_TlOr1))
        print('kb_TlOr2:',"{:.3e}".format(kb_TlOr2))
        print('kb_TrOl1:',"{:.3e}".format(kb_TrOl1))
        print('kb_TrOl2:',"{:.3e}".format(kb_TrOl2))
        print('kb_Td:',"{:.3e}".format(kb_Td))

    rates_backward = {   
            'kb_TR1': kb_TlOr1,
            'kb_TR2': kb_TlOr2,
            'kb_TL1': kb_TrOl1,
            'kb_TL2': kb_TrOl2,
            'kb_Td' : kb_Td
    }
    
    return rates_backward


def get_lida_rates_eyring(Tl,Or1,Or2,breaking_point,invert=False,temperature=TEMPERATURE,verbose=True):
    """
    Here all the kinetic rates are computed as a function of the template
    and oligomer strands given for a LIDA system
    """
    if breaking_point > len(Tl):
        raise BreakPointError

    kf = 2e7
    
    Tr = list(''.join([Or1,Or2]))
    Tl = list(Tl)                       # 5'-3'   NNNNNNNNNNNNNNNNNNNNNN
    Or1 = list(Or1)                     # 3'-5'   NNNNNNNNNN
    Or2 = list(Or2)                     # 3'-5'             NNNNNNNNNNNN
    if invert == True:
        Or1 = list(Or1)[::-1]       
        Or2 = list(Or2)[::-1]           
    Ol1 = Tl[:breaking_point]                          
    Ol2 = Tl[breaking_point:]       

    space1 = ' '.join(['' for i in range(2*len(Or1))])
    space2 = ' '.join(['' for i in range(2*len(Ol1))])
    if verbose == True:
        print(' D: destabilizing \n'
            ' Template:     ',' '.join(Tl),'\n',
            'D-Oligomer1:  ',' '.join(Or1),'\n',
            'D-Oligomer2:  ',space1,' '.join(Or2),'\n',
            'D-Template:   ',' '.join(Tr),'\n',
            'Oligomer1:    ',' '.join(Ol1),'\n',
            'Oligomer2:    ',space2,' '.join(Ol2),'\n',)

    rates = {   
                'forward' :{},
                'backward':{}
            }
    
    #Routines for getting backward reaction rates

    Tl_dx = Tl[:len(Or1)]
    dg_TlOr1 = string_free(Tl_dx,Or1,temperature,verbose=verbose)
    kb_TlOr1 = eyring_backwards(dg_TlOr1,temperature)
    if verbose == True:
        print('TlOr1 Done \n \n')

    Tl_sx = Tl[len(Or1):]
    dg_TlOr2 = string_free(Tl_sx,Or2,temperature,verbose=verbose)
    kb_TlOr2 = eyring_backwards(dg_TlOr2,temperature)
    if verbose == True:
        print('TlOr2 Done\n \n')

    Tr_dx = Tr[:len(Ol1)]
    dg_TrOl1 = string_free(Tr_dx,Ol1,temperature,verbose=verbose)
    kb_TrOl1 = eyring_backwards(dg_TrOl1,temperature)
    if verbose == True:
        print('TrOl1 Done\n \n')

    Tr_sx = Tr[len(Ol1):]
    dg_TrOl2 = string_free(Tr_sx,Ol2,temperature,verbose=verbose)
    kb_TrOl2 = eyring_backwards(dg_TrOl2,temperature)
    if verbose == True:
        print('TrOl2 Done\n \n')

    dg_Td = string_free(Tl,Tr,temperature,verbose=verbose)
    kb_Td = eyring_backwards(dg_Td,temperature)
    if verbose == True:
        print('Td Done\n \n')
    if verbose == True:
        print('kb_TlOr1:',"{:.3e}".format(kb_TlOr1))
        print('kb_TlOr2:',"{:.3e}".format(kb_TlOr2))
        print('kb_TrOl1:',"{:.3e}".format(kb_TrOl1))
        print('kb_TrOl2:',"{:.3e}".format(kb_TrOl2))
        print('kb_Td:',"{:.3e}".format(kb_Td))

    rates_backward = {   
            'kb_TR1': kb_TlOr1,
            'kb_TR2': kb_TlOr2,
            'kb_TL1': kb_TrOl1,
            'kb_TL2': kb_TrOl2,
            'kb_Td' : kb_Td
    }

    return rates_backward


