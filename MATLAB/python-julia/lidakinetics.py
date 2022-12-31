


import values
import math

hs = values.hs
############## Physics Functions

def gibbs_free(enthalpy, entropy, temperature=26): #[DG]=Kcal/mol
    entropy_kcal = entropy/1000 #[S]=cal/mol*K
    return (enthalpy - (273.15+temperature)*entropy_kcal) # [H]=Kcal/mol 
    
def k_equilibrium(DG, temperature=26):
    R = 1.987e-3                                             # [(Kcal)/(mol*k)]
    return math.exp(-(DG/(R*(273.15+temperature))))   # (Kcal/mol)/((Kcal/mol*K)*K) -> adimensional

############## String Reading Functions

#check whether there are mismatches
def iswc(n1,n2):
    if n1 == 'A' and n2 == 'T' or n1 == 'C' and n2 == 'G' or n1 == 'T' and n2 == 'A' or n1 == 'G' and n2 == 'C':
        return True
    else:
        return False

#check whether there are abasic damages

def posab(n1,n2):
    if n1 == 'U':
        return int(0)
    elif n2 == 'U':
        return int(1)
    else:
        return int

def wc(base:str):
    if base == 'A': return str('T')
    if base == 'T': return str('A')
    if base == 'C': return str('G')
    if base == 'G': return str('C')

def modified_matrix(strandA, strandB):
    strands_org = [list(strandA),list(strandB)]
    print('The original strand matrix is \n',strands_org[0],'\n',strands_org[1])
    strands_mod = strands_org
    for i in [0,1]:
        for n in range(len(strands_org[i])):
            if strands_org[i][n] == 'U':
                strands_mod[i][n] = wc(strands_org[i-1][n])
    if strands_mod != strands_org:
        print('The corrected strand matrix is \n',strands_mod[0],'\n',strands_mod[1])
    return strands_mod

def get_flag(strandA, strandB):
    strands = [list(strandA),list(strandB)]
    flag = {'abasic':[],'terminal_mismatches':False,'mismatches':[],'termAT':False}
    for i in [0,1]:
        for n in range(1,len(strands[i])-1):
            if strands[i][n] == 'U':
                flag['abasic'].append([i,n])
        for n in range(len(strands[i])):
            if strands[i][n] != wc(strands[i-1][n]):
                flag['mismatches'].append([i,n])
            if strands[i][0] or strands[i][-1] == 'A' or 'T':
                flag['termAT'] = True
            if strands[i][0] not in [wc(strands[i-1][0]),'U']:
                flag['terminal_mismatches'] = True
            if strands[i][-1] not in [wc(strands[i-1][-1]),'U']:
                flag['terminal_mismatches'] = True
    return flag


def abasic_ddg(context,opposing_nucleotide):
    if context in [['G','U','G'],['T','U','G'],['G','U','T'],['T','U','T']]:
        return values.ab['GUG'][opposing_nucleotide]
    if context in [['C','U','C'],['A','U','C'],['C','U','A'],['A','U','A']]:
        return values.ab['CUC'][opposing_nucleotide]
    else:
        avg = (values.ab['GUG'][opposing_nucleotide]+values.ab['CUC'][opposing_nucleotide])/2
        return avg

def add_initiation(free_energy_list,temperature):
    free_energy_list.append(gibbs_free(hs['initiation']['H'],hs['initiation']['S'],temperature=temperature))
    print('Adding initiation free energy equal to', gibbs_free(hs['initiation']['H'],hs['initiation']['S'],temperature=temperature))


def terminalATpenalty(strands, free_energy_list, temperature):
    terminal = gibbs_free(hs['TerminalATPenalty']['H'],hs['TerminalATPenalty']['S'],temperature=temperature)
    if strands[0][0] == 'A' and strands[1][0] == 'T':
        free_energy_list.append(terminal)
        print('Adding terminal AT penalty for initial bases equal to:',terminal)
    if strands[0][0] == 'T' and strands[1][0] == 'A':
        free_energy_list.append(terminal)
        print('Adding terminal AT penalty for initial bases equal to:',terminal)
    if strands[0][-1] == 'A' and strands[1][-1] == 'T':
        free_energy_list.append(terminal)
        print('Adding terminal AT penalty for final bases equal to:',terminal)
    if strands[0][-1] == 'T' and strands[1][-1] == 'A':
        free_energy_list.append(terminal)
        print('Adding terminal AT penalty for final bases equal to:',terminal)

def get_dangling_shape(strandA,strandB,unionpoint):
    lenA = len(strandA)
    lenB = len(strandB)
    empty = []
    matrix = [list(strandA),list(strandB)]
    # print(matrix)
    lens = [lenA,lenB]
    if lenA == lenB:
        return empty#there's no dangling end if strands have the same length
    else:
        if lenA > lenB:
            l = 0
            s = 1
        if lenA < lenB:
            l = 1
            s = 0
        if unionpoint == -1:
            matrix = [s[::-1] for s in matrix]
            unionpoint = lens[l]-lens[s]
            return [matrix[l][lens[s]],matrix[l][lens[s]-1],matrix[s][-1], '5']
        elif unionpoint == 0:
            return [matrix[l][lens[s]],matrix[l][lens[s]-1],matrix[s][-1], '3']
        elif (unionpoint+lens[s]) == lens[l]:
            return [matrix[l][unionpoint-1],matrix[l][unionpoint],matrix[s][0], '5']
        else:
            return [[matrix[l][unionpoint+lens[s]],matrix[l][unionpoint+lens[s]-1],matrix[s][-1], '3'],  #DX
                    [matrix[l][unionpoint-1],matrix[l][unionpoint],matrix[s][0], '5']]                   #SX
            
def dangling_free_energy(dangling_shape,temperature):
    dang_dg = []
    if len(dangling_shape) == 2:
        for i in range(len(dangling_shape)):
            d = dangling_shape[i]
            HnS = values.dangling[d[-1]][d[0]][d[1]][d[2]]
            gibbs = gibbs_free(HnS['H'],HnS['S'],temperature)
            dang_dg.append(gibbs)
            return sum(dang_dg)

    else:
        d = dangling_shape
        print(d)
        HnS = values.dangling[d[-1]][d[0]][d[1]][d[2]]
        gibbs = gibbs_free(HnS['H'],HnS['S'],temperature)
        dang_dg.append(gibbs)
        return sum(dang_dg)
    

def terminal_mismatch_check(strandA,strandB):
    initial = bool
    final   = bool
    if strandA[0] != wc(strandB[0]):
        initial = True
    if strandA[-1] != wc(strandB[-1]):
        final = True
    return initial, final

def terminal_mismatch_free_energy(strandA,strandB,temperature):
    term_mis_dg = []
    initial, final = terminal_mismatch_check(strandA,strandB)
    if initial: 
        shape1 = [strandA[0],strandA[1],strandB[1],'5']
        shape2 = [strandB[0],strandB[1],strandA[1],'3']
        dg1 = dangling_free_energy(shape1,temperature)
        dg2 = dangling_free_energy(shape2,temperature)
        term_mis_dg.append(dg1+dg2)
    if final:
        shape1 = [strandA[-1],strandA[-2],strandB[-2],'3']
        shape2 = [strandB[-1],strandB[-2],strandA[-2],'5']
        dg1 = dangling_free_energy(shape1,temperature)
        dg2 = dangling_free_energy(shape2,temperature)
        term_mis_dg.append(dg1+dg2)
    return term_mis_dg

def cut_larger(stringA,stringB):
    lenA = len(stringA)
    lenB = len(stringB)
    if lenA == lenB:
        return stringA, stringB
    if lenA > lenB:
        return stringA[:lenB], stringB 
    if lenA < lenB:
        return stringA, stringB[lenB-lenA:]  
##MAIN

def string_free(strandA, strandB, dang=[], temperature=26):
    print('All values are measured in Kcal/mol')
    # print('Dangling end contribution is approximated as an additional AT basepair at the end: https://academic.oup.com/nar/article/28/9/1929/2903888')
    hs = values.hs
    strands = modified_matrix(strandA,strandB)
    original_strands = [list(strandA),list(strandB)]
    flag = get_flag(strandA, strandB)
    print(flag)
    pairs_dg = []
    dangling_dg = []
    term_mism_dg = []
    abasic_penalties = []
    tot_term_mism = 0
    tot_dangling_contribution = 0
    add_initiation(pairs_dg,temperature)
    if flag['termAT'] == True:
        terminalATpenalty(strands, pairs_dg, temperature)
    for i in range(int(len(strands[0]))-1):
        if strands[0][i] == wc(strands[1][i]):
            pairfree = gibbs_free(hs[strands[0][i]][strands[0][i+1]][strands[1][i+1]]['H'],hs[strands[0][i]][strands[0][i+1]][strands[1][i+1]]['S'],temperature)
            pairs_dg.append(pairfree)
            if strands[0][i+1] == wc(strands[1][i+1]):
                print('Delta G for the WC pair',strands[0][i],strands[0][i+1],'-',strands[1][i],strands[1][i+1],':',pairfree)
            else:
                print('Delta G for the mismatch',strands[0][i],strands[0][i+1],'-',strands[1][i],strands[1][i+1],':',pairfree)
        else:
            mismfree = gibbs_free(hs[strands[1][i+1]][strands[1][i]][strands[0][i]]['H'],hs[strands[1][i+1]][strands[1][i]][strands[0][i]]['S'],temperature)
            pairs_dg.append(mismfree)
            print('Delta G for the mismatch',strands[0][i],strands[0][i+1],'-',strands[1][i],strands[1][i+1],':',mismfree)
    if flag['abasic'] != []:
        for ab in flag['abasic']:
            context = [original_strands[ab[0]][ab[1]+j] for j in [-1,0,1]]
            opp_base = original_strands[ab[0]-1][ab[1]]
            abasic_penalties.append(-abasic_ddg(context,opp_base))
            print('Abasic defect found in strand',ab[0],'at position',ab[1])
            print('Context:',''.join(context)+'; Opposing base:',opp_base)
            print('The free energy penalty is',-abasic_ddg(context,opp_base))
    if 'U' in [original_strands[0][0],original_strands[1][0]]:
        abasic_penalties.append(2)
        print('Added a dangling abasic penalty of 2 at initial position',original_strands[0][0] or original_strands[1][0])
    if 'U' in [original_strands[0][-1],original_strands[1][-1]]:
        abasic_penalties.append(2)
        print('Added a dangling abasic penalty of 2 at final position')

    if 'U' not in dang:
        if (dang != []) and (flag['terminal_mismatches'] == False):
                dangling_dg.append(dangling_free_energy(dang,temperature))
                tot_dangling_contribution = sum(dangling_dg)
                print('dangling_contribution = ',dangling_dg)
        if flag['terminal_mismatches'] == True:
            ######### temp fix
            pass
            # term_mism_dg.append(terminal_mismatch_free_energy(strandA,strandB,temperature))
            # tot_term_mism = sum(term_mism_dg)
            # print('Terminal Mismatches contribute is:', tot_term_mism)

    repaired_freenergy = sum(pairs_dg)
    tot_abasic_penalty = sum(abasic_penalties)
    freenergy = repaired_freenergy + tot_abasic_penalty + tot_dangling_contribution + tot_term_mism
    print('Delta G of repaired strand:',repaired_freenergy)
    print('Dangling ends contribution:',tot_dangling_contribution)
    print('Total abasic penalty:',tot_abasic_penalty)
    print('Final Delta G:',freenergy)
    return freenergy

class BreakPointError(Exception):
    def __init__(self, message='Break point cannot be larger than string'):
        self.message = message
        super().__init__(self.message)
    def __str__(self) -> str:
        return f'{self.message}'
    pass

def get_lida_rates(Tl,Or1,Or2,breaking_point,invert=False,temperature=26,kf=2e7):

    if breaking_point > len(Tl):
        raise BreakPointError
   
    Tr = list(''.join([Or1,Or2]))
    Tl = list(Tl)                       # 5'-3'   LLLLLLLLLLLLLLLLLLLLLL
    Or1 = list(Or1)                     # 3'-5'   RRRRRRRRRR
    Or2 = list(Or2)                     # 3'-5'             RRRRRRRRRRRR

    if invert == True:
        Or1 = list(Or1)[::-1]       
        Or2 = list(Or2)[::-1]

                                        # 3'-5'   RRRRRRRRRRRRRRRRRRRRRR
    Ol1 = Tl[:breaking_point]           # 5'-3'   LLLLLLLLLL     
    Ol2 = Tl[breaking_point:]           # 5'-3'             LLLLLLLLLLLL

    space1 = ' '.join(['' for i in range(2*len(Or1))])
    space2 = ' '.join(['' for i in range(2*len(Ol1))])

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

    Tl_sx = Tl[:len(Or1)]
    dang_TlOr1 = get_dangling_shape(Tl,Or1,0)
    print("3' Dangling end shape is:",dang_TlOr1[0],dang_TlOr1[1],'\n',
          '                          ',dang_TlOr1[2])
    Tl_sx_cut, Or1_cut = cut_larger(Tl_sx,Or1)
    dg_TlOr1 = string_free(Tl_sx_cut, Or1_cut,dang_TlOr1,temperature)
    kd_TlOr1 = k_equilibrium(dg_TlOr1)
    kb_TlOr1 = kf/kd_TlOr1
    print('TlOr1 Done \n \n')

    Tl_dx = Tl[len(Or1):]
    dang_TlOr2 = get_dangling_shape(Tl,Or2,-1)
    print('Dangling end shape is:',dang_TlOr2[0],dang_TlOr2[1],'\n',
          '                       ',dang_TlOr2[2])
    Tl_dx_cut, Or2_cut = cut_larger(Tl_dx,Or2)
    dg_TlOr2 = string_free(Tl_dx_cut, Or2_cut,dang_TlOr2,temperature)
    kd_TlOr2 = k_equilibrium(dg_TlOr2)
    kb_TlOr2 = kf/kd_TlOr2
    print('TlOr2 Done\n \n')

    Tr_sx = Tr[:len(Ol1)]
    dang_TrOl1 = get_dangling_shape(Tr[::-1],Ol1[::-1],-1)
    print('Dangling end shape is:',dang_TrOl1[0],dang_TrOl1[1],'\n',
        '                       ',dang_TrOl1[2])
    Tr_sx_cut, Ol1_cut = cut_larger(Tr_sx,Ol1)
    dg_TrOl1 = string_free(Tr_sx_cut, Ol1_cut,dang_TrOl1,temperature=26)
    kd_TrOl1 = k_equilibrium(dg_TrOl1)
    kb_TrOl1 = kf/kd_TrOl1
    print('TrOl1 Done\n \n')

    Tr_dx = Tr[len(Ol1):]
    dang_TrOl2 = get_dangling_shape(Tr[::-1],Ol2[::-1],0)
    print('Dangling end shape is:',dang_TrOl2[0],dang_TrOl2[1],'\n',
        '                       ',dang_TrOl2[2])
    Tr_dx_cut, Ol2_cut = cut_larger(Tr_dx,Ol2)
    dg_TrOl2 = string_free(Tr_dx_cut, Ol2_cut,dang_TrOl2,temperature=26)
    kd_TrOl2 = k_equilibrium(dg_TrOl2)
    kb_TrOl2 = kf/kd_TrOl2
    print('TrOl2 Done\n \n')

    dang_Orl1 = get_dangling_shape(Or1,Ol1,0)
    if dang_Orl1 != []:
        print('Dangling end shape is:',dang_Orl1[0],dang_Orl1[1],'\n',
              '                       ',dang_Orl1[2])
    if dang_Orl1 == []:
        print('Oligomers 1 are of the same length')
    Or1_cut, Ol1_cut = cut_larger(Or1,Ol1)
    dg_Orl1 = string_free(Or1_cut,Ol1_cut,dang_Orl1,temperature=26)
    kd_Orl1 = k_equilibrium(dg_Orl1)
    kb_Orl1 = kf/kd_Orl1
    print('Orl1 Done\n \n')

    dang_Orl2 = get_dangling_shape(Or2,Ol2,-1)
    if dang_Orl2 != []:
        print('Dangling end shape is:',dang_Orl2[0],dang_Orl2[1],'\n',
            '                       ',dang_Orl2[2])
    if dang_Orl2 == []:
        print('Oligomers 1 are of the same length')
    Or2_cut, Ol2_cut = cut_larger(Or2,Ol2)
    dg_Orl2 = string_free(Or2_cut,Ol2_cut,dang_Orl2,temperature=26)
    kd_Orl2 = k_equilibrium(dg_Orl2)
    kb_Orl2 = kf/kd_Orl2
    print('Orl2 Done\n \n')
    

    dg_Td = string_free(Tl,Tr,temperature=26)
    kd_Td = k_equilibrium(dg_Td)
    kb_Td = kf/kd_Td
    print('Td Done\n \n')

    print('kb_TlOr1:',"{:.3e}".format(kb_TlOr1))
    print('kb_TlOr2:',"{:.3e}".format(kb_TlOr2))
    print('kb_TrOl1:',"{:.3e}".format(kb_TrOl1))
    print('kb_TrOl2:',"{:.3e}".format(kb_TrOl2))
    print('kb_Orl1:',"{:.3e}".format(kb_Orl1))
    print('kb_Orl2:',"{:.3e}".format(kb_Orl2))
   
    print('kb_Td:',"{:.3e}".format(kb_Td))

    DHK = { 'kb_TR1'    :kb_TlOr1,
            'kb_TR2'    :kb_TlOr2,
            'kb_TL1'    :kb_TrOl1,
            'kb_TL2'    :kb_TrOl2,
            'kb_O1'     :kb_Orl1,
            'kb_O2'     :kb_Orl2,
            'kb_duplex' :kb_Td}

    EQK = { 'ke_TR1'    :kd_TlOr1,
            'ke_TR2'    :kd_TlOr2,
            'ke_TL1'    :kd_TrOl1,
            'ke_TL2'    :kd_TrOl2,
            'ke_O1'     :kd_Orl1,
            'ke_O2'     :kd_Orl2,
            'ke_duplex' :kd_Td}
    
    return DHK, EQK

def get_vanilla_rates(str1,str2,temperature=26,kf=2e7):
    cut1, cut2 = cut_larger(str1,str2)
    dg = string_free(cut1,cut2,dang=[],temperature=temperature)
    kd = k_equilibrium(dg)
    kb = kf/kd
    return kd, kb



#############################  Deprecated stuff

        # if strands[0][i+1] == 'U':
        #         context = [strands[0][j] for j in [i,i+1,i+2]]
        #         pairs_dg.append(abasic_ddg(context,strands[1][i+1]))
        #         print('Abasic defect in strand 1 at position:', i+1)

        # if strands[1][i+1] == 'U':
        #         context = [strands[1][j] for j in [i,i+1,i+2]]
        #         pairs_dg.append(abasic_ddg(context,strands[0][i+1]))
        #         print('Abasic defect in strand 2 at position:', i+1)


        # if strands[0][i+1] or strands[1][i+1] == 'U':
        #         num = posab(strands[0][i+1],strands[1][i+1])
        #         context = [strands[num][j] for j in [i,i+1,i+2]]
        #         pairs_dg.append(abasic_ddg(context,strands[num-1][i+1]))
        #         skip = 1
        #         continue