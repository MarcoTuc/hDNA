import networkx as nx
import copy
import pandas as pd 

from classes.chamber import Chamber
from classes.model import Model
from classes.strand import Strand
 
class Kinetwork(object):

    """Given a chamber it will generate the corresponding
    Kinetic network to be later passed to the simulator"""

    def __init__(self, model: Model, s1: Strand, s2: Strand, minimum_nucleation):

        self.model = model 
        self.s1 = s1
        self.s2 = s2
        self.mincore = minimum_nucleation
        self.chamber = Chamber(self.model, self.s1, self.s2, self.mincore)


        self.Graph = nx.Graph()
        self.add_nodes()
        self.add_reactions()
        self.clean_duplicates()


    def add_nodes(self):

        self.Graph.add_node(self.chamber.ssstruct, structure = self.chamber.ssstruct, state = 'singlestranded', pairs = 0)

        for s in self.chamber.offcores:
            self.Graph.add_node(s, object = s, structure = s.structure, state = 'off_register', pairs = int(s.total_nucleations))

        for s in self.chamber.oncores:
            self.Graph.add_node(s, object = s, structure = s.structure, state = 'on_register', pairs = int(s.total_nucleations))

        for s in self.chamber.oncores:
            for z in s.zipping[1:]: #first zipping value is the oncore itself so flush it away by indexing
                self.Graph.add_node(z, structure = z.structure, state = 'zipping', pairs = int(z.total_nucleations))

        self.Graph.add_node(self.chamber.duplex, structure = self.chamber.duplex.structure, state = 'duplex', pairs = self.chamber.duplex.total_nucleations)
        """   note that duplicated nodes will not be added to the resulting graph, hence I can run through 
              all zipping states and be sure that I will only get one note per each zipping that is common
              from different starting native nucleations. """

    def add_reactions(self):

        ON  = self.node_filter('state','on_register')
        OFF = self.node_filter('state','off_register')
        D = self.node_filter('state', 'duplex')
        print(D)

        for on, data in ON.nodes.items():
            self.Graph.add_edge(self.chamber.ssstruct, on, kind = 'on_nucleation') #ADDPROPERTY
            self.Graph.add_edge(on, on.zipping[1], kind = 'zipping')
            for i, z in enumerate(on.zipping[2:], start=2):
                self.Graph.add_node(z, structure = z.structure, state = 'zipping', pairs = int(z.total_nucleations))
                self.Graph.add_node(on.zipping[i-1], structure = on.zipping[i-1].structure, state = 'zipping', pairs = int(z.total_nucleations))
                self.Graph.add_edge(z, on.zipping[i-1], kind = 'zipping')
            self.Graph.add_edge(on.zipping[-1], self.chamber.duplex, kind = 'zipping-end')

        L, R = self.chamber.split_offcores()

        for i, (left, right) in enumerate(zip(L, R)):
            self.Graph.add_edge(self.chamber.ssstruct, left, kind = 'off_nucleation')
            self.Graph.add_edge(self.chamber.ssstruct, right, kind = 'off_nucleation')
            if i > 0: 
                self.Graph.add_edge(left, L[i-1], kind = 'sliding')
                self.Graph.add_edge(right, R[i-1], kind = 'sliding')
        self.Graph.add_edge(L[-1], list(D.nodes())[0], kind = 'sliding-end')
        self.Graph.add_edge(R[-1], list(D.nodes())[0], kind = 'sliding-end')

    def clean_duplicates(self):
        """ Unfortunately my code makes duplicates.
            Also unfortunately it was easier to just get rid of them
            after they are made than fix the entire code """
        df = pd.DataFrame(self.Graph.nodes(data=True))
        labels = df[0]
        struct = [list(df[1])[i]['structure'] for i in range(len(list(df[1])))]
        diz = dict(zip(labels, struct))
        Relabeled = nx.relabel_nodes(self.Graph, diz)
        self.Graph = Relabeled

        #TODO: FIX THIS SHIT DOWN HERE  

        # for off, data in OFF.nodes.items():
        #     self.Graph.add_edge(self.chamber.ssstruct, off, kind = 'off_nucleation') #ADDPROPERTY
        #     offleft, offright = self.chamber.split_offcores()
        #     for i, (l, r) in enumerate(zip(offleft, offright)):
        #         if 0 < i < len(offleft): 
        #             self.Graph.add_edge(l, offleft[i-1].structure, kind = 'sliding') #ADDPROPERTY
        #             self.Graph.add_edge(r, offright[i-1].structure, kind = 'sliding') #ADDPROPERTY
        #     #   TODO
        #     #   add here a condition for sliding into duplex from the most duplexed sliding state 
        #     self.Graph.add_edge(offleft[-1], self.chamber.duplex, kind = 'sliding_closure') #ADDPROPERTY
        #     self.Graph.add_edge(offright[-1], self.chamber.duplex, kind = 'sliding_closure') #ADDPROPERTY
    
    
    def node_filter(self, property, attribute):
        return self.Graph.subgraph( 
        [n for n, attrdict in self.Graph.nodes.items() if attrdict [str(property)] == str(attribute)])

    def get_neighbor_zippings(self, structure, onlyup = False, onlydown = False):

        left, right = structure.split('+')

        def get_i(lst):
            indices = []
            for i, el in enumerate(lst):
                if i > 0 and el != lst[i-1]:
                    indices.append(i)
            return indices

        def update_structure(string, character: str):
            indices = get_i(string)
            indices_inv = get_i(string[::-1])
            updated = string
            for index in indices:
                updated = updated[:index-1] + character + updated[index:]
            updated_inv = updated[::-1]
            for index in indices_inv:
                updated_inv = updated_inv[:index-1] + character + updated_inv[index:]
            return updated_inv[::-1]

        up = '+'.join([update_structure(left,'('),update_structure(right,')')])
        down = '+'.join([update_structure(left,'.'),update_structure(right,'.')])

        if onlyup == True:
            return up
        elif onlydown == True:
            return down
        else: return up, down

# DEPRECATED (?)

    # def make_sliding_transitions(self):

    #     left_off = []
    #     right_off = []
    #     for off in range(len(self.offnodes)):
    #         if self.offnodes[off].offregister == 'left':  left_off.append(self.offnodes[off]) 
    #         if self.offnodes[off].offregister == 'right': right_off.append(self.offnodes[off]) 
        
    #     left_sliding_transitions = []
    #     left_sliding_nucleations = []
    #     for sliding in left_off:
            
    #         left_sliding_nucleations.append
    #         #TODO this dictionary of reaction rules 
    #         (
    #             {
    #             'forward'   :{
    #                     'name':None,
    #                     'rate':None,
    #                     'rule':None},
    #             'backward'  :{
    #                     'name':None,
    #                     'rate':None,
    #                     'rule':None}
    #             }
    #         )
    #         left_sliding_transitions.append
    #         #TODO this dictionary of reaction rules 
    #         (
    #             {
    #             'forward'   :{
    #                     'name':None,
    #                     'rate':None,
    #                     'rule':None},
    #             'backward'  :{
    #                     'name':None,
    #                     'rate':None,
    #                     'rule':None}
    #             }
    #         )

    #     right_sliding_nucleations = []
    #     right_sliding_transitions = []
    #     for sliding in right_off:

    #         right_sliding_nucleations.append
    #         #TODO this dictionary of reaction rules 
    #         (
    #             {
    #             'forward'   :{
    #                     'name':None, #put here singlestrands --> first_offregister_nucleation_structure AS A NAME 
    #                     'rate':None, #put here kf_nucleation_geometric as a general parameter (later I may make this a function of the number of nucleotides of separation)
    #                     'rule':None},#put here SS + SS --> offregister_nucleation_structure
    #             'backward'  :{
    #                     'name':None, #put here first_offregister_nucleation_structure --> singlestrands
    #                     'rate':None, #put here kf_nucleation_geometric/K_eq (compute K_eq by taking the free energy difference between origin and destination structure and plugging it into a K_eq function)
    #                     'rule':None} #put here offregister_nucleation_structure --> SS + SS
    #             }
    #         )

    #         right_sliding_transitions.append
    #         #TODO this dictionary of reaction rules 
    #         (
    #             {
    #             'forward'   :{
    #                     'name':None, #put here struct_origin --> struct_destin AS A NAME 
    #                     'rate':None, #put here kf_sliding as a general parameter (later I may make this a function of the number of nucleotides of separation)
    #                     'rule':None},#put here struct_origin --> struct_destin
    #             'backward'  :{
    #                     'name':None, #put here struct_destin --> struct_origin
    #                     'rate':None, #put here kf_sliding/K_eq (compute K_eq by taking the free energy difference between origin and destination structure and plugging it into a K_eq function)
    #                     'rule':None} #put here struct_destin --> struct_origin AS A REACTION RULE 
    #             }
    #         )

    # def make_zipping_transitions(self):
    #     native_nucleations = []
    #     zipping_transitions = []
    #     for on in self.onnodes:
    #         #PSEUDO: append native nucleations 
    #         for z in on.zipping:
    #             zipping_transitions.append
    #             (
    #                 {
    #                 'forward'   :{
    #                         'name':None, #put here struct_origin --> struct_destin AS A NAME 
    #                         'rate':None, #put here kf_zipping as a general parameter (later I may make this a function of the number of nucleotides of separation)
    #                         'rule':None},#put here struct_origin --> struct_destin
    #                 'backward'  :{
    #                         'name':None, #put here struct_destin --> struct_origin
    #                         'rate':None, #put here kf_zipping/K_eq (compute K_eq by taking the free energy difference between origin and destination structure and plugging it into a K_eq function)
    #                         'rule':None} #put here struct_destin --> struct_origin AS A REACTION RULE 
    #                 }
    #             )

        