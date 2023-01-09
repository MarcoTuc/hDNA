import networkx as nx

from classes.chamber import Chamber
from classes.model import Model
 
class Kinetwork(object):

    """Given a chamber it will generate the corresponding
    Kinetic network to be later passed to the simulator"""

    def __init__(self, model: Model, chamber: Chamber, minimum_nucleation):

        self.model = model 
        self.chamber = chamber
        self.mincore = minimum_nucleation

        self.Graph = nx.Graph()
        self.add_nodes()

    def add_nodes(self):
        self.Graph.add_node(self.chamber.ssstruct, state = 'singlestranded')

        for s in self.chamber.offcores:
            self.Graph.add_node(s.structure, state = 'off_register', pairs = s.total_nucleations)
        for s in self.chamber.oncores:
            self.Graph.add_node(s.structure,  state = 'on_register', pairs = s.total_nucleations)
        for z in self.chamber.oncores:
            for s in z.zipping: 
                self.Graph.add_node(s.structure, state = 'zipping', pairs = s.total_nucleations)

        self.Graph.add_node(self.chamber.duplex.structure, state = 'duplex', pairs = self.chamber.duplex.total_nucleations)
        # note that duplicated nodes will not be added to the resulting graph, hence I can run through 
        # all zipping states and be sure that I will only get one note per each zipping that is common
        # from different starting native nucleations. 

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

        