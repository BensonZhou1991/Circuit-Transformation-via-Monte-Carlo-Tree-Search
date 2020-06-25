# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 19:19:04 2019

@author: zxz58
"""
import sys
sys.path.append('..')
import networkx as nx
from networkx import DiGraph
from operation import FindExecutableNode, ExecuteAllPossibileNodesInDG
from operation import FindInvolvedPhyQInOneLayer
from operation import SWAPInArchitectureGraph
from operation import GetSubCxList
import numpy as np
from qiskit import QuantumCircuit
from qiskit import QuantumRegister
import time
import matplotlib.pyplot as plt

bias_for_score = 0 #default 0
score_decay_rate = 0.7 #default 0.7 for IBM Q20, J-P
visit_decay_rate = 50 # recommended value is 10 or 50?
visit_punish = 0 # whether we use punish visit time brought by back propagation
SWAP_cost = 3
display_state = 1 #whether we print the internal states during search process
delete_trivival_swap = 0 #should we ignore duplicated swaps to avoid dead loop?
delete_trivival_swap_decision = 1 # recommend to be 1 to avoid dead loop
fallback_value = 15
log_data = 0 # transcribe relevant data for figuring?
# deactivation control
max_depth = 2000 # we don't wanna selection procedure to explore too deeply
max_visit = 50000 # we don't wanna same node to be selected too many times
min_deactivation_score = 0.5 # default 0.5
use_heuristic_score = 0

'''
select_mode specifies the mode for evaluation during selection,
['our method'] is a randomly selecting method
['KS', arg] is a deterministic method
['KS_random', (arg, power_for_random)] is a random method 
['ZHOU', (10, -10)]
['ZHOU2', (T_h, T_v)]
['ZHOU3', (bias, T_v)]
['cx_vs_visit', arg]
recommended value:
    Q20: ['KS', 20] or ['our method']
    J_P: ['KS', 5] or ['ZHOU', (5, 1/80 for small circuit or 1/20 for large)]
    A-B-S: ['KS', 5]
'''
#_select_mode_ = ['ZHOU', (5, 1/20)]
_select_mode_ = ['KS', 20]
#_select_mode_ = ['KS_random', (5, 5)]
#_select_mode_ = ['ZHOU3', (6, 1/20)]

'''
mode for Back Propagation
['globalscore_depth', [score_decay_rate, depth_threshold]]
recommended value:
    Q20: ['globalscore', [score_decay_rate]]
    J_P: ['globalscore', [score_decay_rate]]
'''
mode_BP = ['globalscore', [score_decay_rate]]
#mode_BP = ['globalscore_depth', [score_decay_rate]]
#mode_BP = ['globalscore_and_ave_cx', [2.5]]#2.5
#mode_BP = ['globalscore_modified', [score_decay_rate]]

'''
mode for decision
    ['global_score']:
        only consider global score
    ['sim_fix_CX', decay_score, decay_sim]:
        consider both global score and swaps got by simulation
        decay is the weight for simulation scores
        recommended value: ['sim_fix_CX', 1, 0.7] or ['sim_fix_CX', 1, 1.2]
    ['num_sim_swap']
'''
#mode_decision = ['num_sim_swap']
mode_decision = ['global_score']
#mode_decision = ['sim_fix_CX', 1, 1.2]

'''
WARNING: if you add more than one swap in a new node, the following function
will execute all these swaps at one time and renew mapping and vertex list, it
may be not optimal 
'''
def GetNewMapFromFather(current_node, MCTree):
    '''get new mapping via from father node'''
    father_node = MCTree.nodes[current_node]['father_node']
    mapping = MCTree.nodes[father_node]['mapping'].Copy()
    return mapping

def GetNewMap(current_node, MCTree, swap):
    '''update new mapping via exe single swap'''
    mapping = MCTree.nodes[current_node]['mapping']
    mapping.RenewMapViaExchangeCod(swap[0], swap[1])
    return mapping

def RenewExecutionList(current_node, MCTree):
    '''
    execute CX and update lists executed_vertex and executable_vertex
    according to current mapping of a given node (current_node)
    Return:
        number of newly executed CNOT gates
    '''
    DG = MCTree.DG
    AG = MCTree.AG
    father_node = MCTree.nodes[current_node]['father_node']
    if MCTree.nodes[current_node]['executed_vertex'] == None:
        '''this is the first swap among all swaps to be added'''
        if father_node != None: 
            executed_vertex = MCTree.nodes[father_node]['executed_vertex'].copy()
            executable_vertex = MCTree.nodes[father_node]['executable_vertex'].copy()
        else:
            executed_vertex = []
            executable_vertex = FindExecutableNode(DG)
    else:
        executed_vertex = MCTree.nodes[current_node]['executed_vertex']
        executable_vertex = MCTree.nodes[current_node]['executable_vertex']
    num_executed_CX_before = len(executed_vertex)
    mapping = MCTree.nodes[current_node]['mapping']
    res = ExecuteAllPossibileNodesInDG(executable_vertex, executed_vertex, AG, 
                                       DG, mapping,
                                       out_removed_nodes=True,
                                       draw=False, DiG=None, edges_DiG=None)
    executed_vertex, executable_vertex, executed_vertex_current = res
    MCTree.nodes[current_node]['executed_vertex'] = executed_vertex
    MCTree.nodes[current_node]['executable_vertex'] = executable_vertex
    MCTree.nodes[current_node]['executed_vertex_current'] = executed_vertex_current
    num_executed_CX_after = len(executed_vertex)
    return num_executed_CX_after - num_executed_CX_before
    


def CalScore(MCTree, node):
    '''Calculate local score and store it in the node'''
    #cost_g = MCTree.nodes[node]['num_total_add_gates']
    #cost_h = MCTree.nodes[node]['cost_h']
    num_g = len(MCTree.nodes[node]['executed_vertex'])
    #num_h = MCTree.nodes[node]['num_remain_vertex']
    '''different ways for calculating score'''
    #score = (num_g + num_h)/(cost_h + cost_g) #use both existing and heuristic cost
    #score = (num_g) / (cost_g + 0.001) #use only existing cost
    father_node = MCTree.nodes[node]['father_node']
    score = num_g - len(MCTree.nodes[father_node]['executed_vertex'])#executed CNOT only
    return score
    
def ConductCNOTInDGAlongPath(DG, vertex, path, mapping,
                             executable_vertex=None,
                             executed_vertex=None, G=None):
    '''
    conduct CNOT in a vertex in DG along a specific path([control, ..., target])
    in architecture graph, then, renew physical circuit and mapping
    input:
        remove_node: will the node in DG be removed?
        exe_other_CX: should we check and execute other possible CNOT whem swapping?
                      ATTENTION: this fuction is only valid for undirected AG
    '''
    flag_CX = []
    add_gates_count = 0
    v_c_pos = 0
    v_t_pos = len(path) - 1
    num_swaps = len(path) - 2
    flag_head = True #decide in which way (control to target or target to control) we perform SWAP

    for i in range(num_swaps):
        add_gates_count += 3
        if flag_head == True:
            SWAPInArchitectureGraph(path[v_c_pos], path[v_c_pos+1], mapping,
                                    q_phy=None, cir_phy=None, draw=False)
            v_c_pos += 1
            flag_head = not flag_head
        else:
            SWAPInArchitectureGraph(path[v_t_pos], path[v_t_pos-1], mapping,
                                    q_phy=None, cir_phy=None, draw=False)
            v_t_pos -= 1
            flag_head = not flag_head
        res = ExecuteAllPossibileNodesInDG(executable_vertex, executed_vertex,
                                           G, DG, mapping, draw=False, DiG=None,
                                           edges_DiG=None, cir_phy=None,
                                           q_phy=None, out_removed_nodes=True)
        executed_vertex, executable_vertex, removed_nodes = res
        flag_CX.append(len(removed_nodes))
        
    return executable_vertex, executed_vertex, flag_CX

def FindSWAPsInvolved(swaps, executable_vertex, DG, mapping):
    '''
    return swaps having effect on any executable gate
    '''
    q_phy = []
    inv_swaps = []
    for vertex in executable_vertex:
        op = DG.nodes[vertex]['operation']
        q0 = op.involve_qubits[0]
        q1 = op.involve_qubits[1]
        v0 = mapping.DomToCod(q0)
        v1 = mapping.DomToCod(q1)
        q_phy.extend([v0, v1])
    for swap in swaps:
        swap1 = swap[0]
        if (swap1[0] in q_phy) or (swap1[1] in q_phy):
            inv_swaps.append(swap)
    return inv_swaps

class MCTree(DiGraph):
    def __init__(self, AG, DG, shortest_length_AG, shortest_path_AG,
                 swap_combination, mode_sim):
        '''
        swap_combination is a list of swaps to be considered
        T: ratio for node evaluation
        node_count: index for newly added node
        node_size: size for the tree, this variable should be updated continuously
        '''
        super().__init__()
        self.p_swaps = []
        self.node_count = 0
        self.node_size = 0
        self.best_node = None
        self.swap_combination = swap_combination
        self.AG = AG
        self.num_q_log = 16
        self.DG = DG
#        self.search_depth = np.zeros(50) #record depth for last 50 selections
#        self.depth_pos = -1
        self.ave_depth = 0
        self.min_depth = 0
        self.num_CX = len(DG.nodes())
        self.shortest_length_AG = shortest_length_AG
        self.shortest_path_AG = shortest_path_AG
        #self.max_shortest_length_AG = max(shortest_length_AG)
        self.best_score_total = None
        self.root_node = 0
        self.T = visit_decay_rate
        self.start_time = time.time()
        self.first_finish_node = None
        self.first_finish_time = None
        self.first_finish_add_swaps = None
        self.end_time = None
        self.mode_sim = mode_sim
        self.fallback_value = fallback_value
        self.fallback_count = 0
        self.select_mode = _select_mode_
        self.deactivate_nodes = []
        self.simulated_nodes = [] # store simulated nodes in current decision process
        self.log = {'score': [],
                    'visit': [],
                    'swap_vs_cx': [],
                    'num_remain_cx': []} # used to store some data for future analysis
    
    def AddNode(self, father_node, added_SWAP, map_new=None):
        '''
        visited_time: how many times this node has been visited
        best_score: best score for all simulations
        mapping: current mapping
        num_total_add_gates: total number of added (swap) gates
        added_SWAP: list like [(in1, in2),...]
        executed_vertex_current: newly executed nodes in DG
        executed_vertex: all executed nodes in DG
        executable_vertex: nodes in DG for the front layer
        heuristic_score: heuristic score for added node
        heuristic_cost: heuristic cost for added node
        son_nodes: list of indexed of son nodes
        swap_decay: deacy parameter for added SWAP, this is only uesful for 
                    modified globalscore BP mode
        num_remain_vertex: number of unexecuted CNOT gates in logical circuit
        score: local score
        global_score: global score considering all its son nodes, the initial
                      value is its (local) score
        visited_time_total: for some backpropagate protocol, we only add visit_time
                            conditionallyh, hence, this varible is to record the
                            real visited time
        
        return: generated node number
        '''
        new_node = self.node_count
        self.node_count += 1
        self.node_size += 1
        self.add_node(new_node)
        self.nodes[new_node]['father_node'] = father_node
        self.nodes[new_node]['son_nodes'] = []
        self.nodes[new_node]['executed_vertex'] = None
        self.nodes[new_node]['sim_score'] = 0
        self.nodes[new_node]['score'] = 0 #+ 0.05 # We add 0.1 as a bias to improve performance
        self.nodes[new_node]['heuristic_score'] = None
        self.nodes[new_node]['heuristic_cost'] = 0
        self.nodes[new_node]['swap_decay'] = None
        self.nodes[new_node]['activation'] = 1
        self.nodes[new_node]['max_weight'] = 1 # for cal h score
        self.nodes[new_node]['search_depth'] = []
        if father_node != None:
            self.add_edge(father_node, new_node)
            map_new = GetNewMapFromFather(new_node, self)
            self.nodes[new_node]['mapping'] = map_new
            self.nodes[new_node]['num_total_add_gates'] = (self.nodes[father_node]['num_total_add_gates']
                                                            + len(added_SWAP))
        else:
            '''added node is the initial (root) node'''
            self.nodes[new_node]['mapping'] = map_new
            self.nodes[new_node]['num_total_add_gates'] = 0
            add_score = RenewExecutionList(new_node, self)
            self.nodes[new_node]['score'] = add_score
            
        '''update mapping, execut lists and scores'''
        decay = 1
        for swap in added_SWAP:
            map_new = GetNewMap(new_node, self, swap)
            self.nodes[new_node]['mapping'] = map_new
            add_score = RenewExecutionList(new_node, self)
            '''update (local) score for added nodes'''
            self.nodes[new_node]['score'] += add_score * decay
            decay = decay * score_decay_rate
        '''initialize global score'''
        self.nodes[new_node]['global_score'] = self.nodes[new_node]['score']
        
        self.nodes[new_node]['added_SWAP'] = added_SWAP
        self.nodes[new_node]['visited_time'] = 0
        self.nodes[new_node]['visited_time_total'] = 0
        #self.nodes[new_node]['rewards'] = 0
        #self.nodes[new_node]['global_score'] = None
        self.nodes[new_node]['num_remain_vertex'] = self.num_CX - len(self.nodes[new_node]['executed_vertex'])
        if self.nodes[new_node]['num_remain_vertex'] == 0:
            '''if it is finished, then renew relevant variables'''
#            print(self.first_finish_time)
            if self.first_finish_node == None:
                self.first_finish_node = new_node
                self.first_finish_time = time.time()
                self.first_finish_add_swaps = self.nodes[new_node]['num_total_add_gates']
        
        return new_node

    def UpdateSearchDpeth(self, new_depth, fist_sele_node):
        '''add serach depth of each selection to son of root node'''
# =============================================================================
#         pos = self.depth_pos
#         next_pos = pos + 1
#         if next_pos > 49: next_pos = 0
#         self.depth_pos = next_pos
# =============================================================================
        self.last_depth = new_depth
        self.nodes[fist_sele_node]['search_depth'].append(new_depth)
        
    def GetLastSearchDpeth(self):
        #return self.search_depth[self.depth_pos]
        return self.last_depth
    
    def CalWeightedMatrix(self):
        '''
        Calculate weights for each logical qubits pair, and add the single
        weight to eahc CNOT gate operation
        'method' for InitialCostMatrixWeighted
            2: lineal decay
            4: No decay
        '''
        from inputs.map import InitialCostMatrixWeighted
        cost_m = InitialCostMatrixWeighted(self.DG.copy(),
                                           self.num_q_log,
                                           method=2,
                                           add_weight=True)
        '''transform to upper matrix'''
        for raw in range(self.num_q_log):
            for col in range(raw):
                cost_m[col][raw] += cost_m[raw][col]
                cost_m[raw][col] = 0
        self.nodes[0]['weighted_cost_matrix'] = cost_m
        return cost_m
    
    def DeactivateNode(self, node):
        self.nodes[node]['activation'] = 0
        self.deactivate_nodes.append(node)
        '''check whether all son nodes are deactivated'''
        father = self.nodes[node]['father_node']
        _, activations = self.GetSonAttributes(father, ['activation'])
        if np.sum(activations) == 0 and node != self.root_node:
            '''if yes, we deactivate its father node'''
            self.DeactivateNode(father)

    def ResetActivation(self, mode=['all']):
        if mode[0] == 'all':
            while len(self.deactivate_nodes) != 0:
                node = self.deactivate_nodes.pop()
                self.nodes[node]['activation'] = 1
        if mode[0] == 'son of root':
            sons = self.nodes[self.root_node]['son_nodes']
            for son in sons:
                if self.nodes[son]['activation'] == 0:
                    self.nodes[son]['activation'] = 1
                    self.deactivate_nodes.remove(son)
        if mode[0] == 'son':
            sons = self.nodes[mode[1]]['son_nodes']
            for son in sons:
                if self.nodes[son]['activation'] == 0:
                    self.nodes[son]['activation'] = 1
                    self.deactivate_nodes.remove(son)
    
    def CalScore(self, node):
        '''
        Calculate initial score for new node, notice that this score may be 
        renewed during search process
        '''
        print('discarded function!!!')
        score = CalScore(self, node)
        self.nodes[node]['score'] += score
        self.nodes[node]['global_score'] += score
        return score
    
    def AddHScoreToScore(self, node, arg):
        '''
        add h score to the local or global score (represented by arg) of node
        arg: 'score' or 'global_score'
        '''
        print('discarded function!')
        h_score = self.GetHeuristicScore(node)
        h_score = h_score/3
        #print('added h score is', h_score)
        #print('added before is', self.nodes[node]['score'])
        self.nodes[node][arg] += h_score
        #print('added after is', self.nodes[node]['score'])
    
    def FindNonTrivivalSWAPs(self, node):
        swaps = []
        executable_vertex_current = self.nodes[node]['executable_vertex']
        current_map = self.nodes[node]['mapping']
        # New method
        swaps = FindSWAPsInvolved(self.swap_combination,
                                  executable_vertex_current,
                                  self.DG, current_map)        
# =============================================================================
#         # Old slow method
#         for swap in self.swap_combination:
#             if len(swap) == 1:
#                 single_swap = swap
#             else:
#                 '''
#                 if this swap contains multi swaps, we use the first and last
#                 swap to judge whether it is trivival
#                 '''
#                 single_swap = [(swap[0][0], swap[1][0])]
#             flag_nontrivial = CheckSWAPInvolved(single_swap,
#                                                 executable_vertex_current,
#                                                 self.DG, current_map)
#             if flag_nontrivial == True: swaps.append(swap)
# =============================================================================
        return swaps
    
    def FindNonTrivivalSWAPsLi(self, node):
        '''
        this version is based on Pro Li's suggestion
        only support single swap
        '''
        swaps = []
        executable_vertex_current = self.nodes[node]['executable_vertex'].copy()
        #add nodes in the next level
        for vertex in executable_vertex_current:
            next_gates = self.DG.succ[vertex]
            for vertex_next in next_gates:
                if not vertex_next in executable_vertex_current:
                    executable_vertex_current.append(vertex_next)
        
        current_map = self.nodes[node]['mapping']
        swaps = FindSWAPsInvolved(self.swap_combination,
                                  executable_vertex_current,
                                  self.DG, current_map)
        
        return swaps
    
    def ExpandNodeViaSWAP(self, node, swap):
        added_node = self.AddNode(node, swap)
        '''cal h score, it can be commented if necessary'''
        if _select_mode_[0] == 'ZHOU':
            self.AddHScoreToScore(added_node, 'score')
        if _select_mode_[0] == 'ZHOU2' or _select_mode_[0] == 'ZHOU3':
            self.GetHeuristicScore2(added_node)
        if _select_mode_[0] == 'KS' or _select_mode_[0] == 'KS_random': 
            if mode_BP[0] == 'globalscore_and_ave_cx':
                self.GetHeuristicScoreAveSWAP(added_node, mode_BP[1][0])
            else:
                if use_heuristic_score == True:
                    self.GetHeuristicScore3(added_node)
        '''add new node number to its father node's son node list'''
        self.nodes[node]['son_nodes'].append(added_node)
        '''calculate heuristic cost and score for added nodes'''
        # we comment it because we have done this in AddNode function
# =============================================================================
#         #self.CalHeuristicCost(added_node)
#         self.CalScore(added_node)
# =============================================================================
        return added_node
    
    def ExpandNodeViaSWAPs(self, node, SWAPs):
        '''expand a node via SWAP list, [swap1, swap2, ...]'''
        added_nodes = []
        '''the following for loop may be improved via map'''
        for swap in SWAPs:
            added_nodes.append(self.ExpandNodeViaSWAP(node, swap))
        return added_nodes
    
    def ExpandNode(self, node):
        '''
        expand a node via all non-trivival swaps
        '''
        if self.out_degree[node] != 0: raise(Exception('Expanded node already has son nodes.'))
        swaps = self.FindNonTrivivalSWAPs(node)
        #swaps = self.FindNonTrivivalSWAPsLi(node)
        self.p_swaps.append(len(swaps))
        
        if delete_trivival_swap == 1:
            '''
            check wehther expanded node does not exe any CX gates, if yes, we delete
            the swap that the father of expanded node does to avoid dead loop
            '''
            if node != 0:
                if self.nodes[node]['score'] == 0:
                    swap_delete =  self.nodes[node]['added_SWAP']
                    if swap_delete in swaps: swaps.remove(swap_delete)
        
        added_nodes = self.ExpandNodeViaSWAPs(node, swaps)
        return added_nodes
    
    def GetSonScores(self, node, arg='score'):
        sons = self.nodes[node]['son_nodes']
        scores = np.empty([len(sons)])
        visited_time = np.empty([len(sons)])
        pos = -1
        for son in sons:
            pos += 1
            scores[pos] = (self.nodes[son][arg])
            visited_time[pos] = (self.nodes[son]['visited_time'])
        return scores, visited_time, sons

    def GetSonAttributes(self, node, args):
        '''get attributes and sons, represented in list args, from all sons of node'''
        sons = self.nodes[node]['son_nodes']
        num_attr = len(args)
        res = []
        for _ in range(num_attr): res.append(np.empty([len(sons)]))
        pos_son = -1
        for son in sons:
            pos_son += 1
            for pos_arg in range(num_attr):
                new_value = self.nodes[son][args[pos_arg]]
                if new_value == None:
                    #self.PrintSonNodesArgs(node, args)
                    #raise(Exception('Value None'))
                    new_value = 0;
                res[pos_arg][pos_son] = new_value
        return sons, res
    
    def PickBestSonNodeRandomly(self, node, arg):
        '''this is a subfunction for selection section'''
        select_mode = self.select_mode
        if select_mode[0] == 'our method':
            '''arg is the criterion argument'''
            scores, visited_time, sons = self.GetSonScores(node, arg)
            scores = np.exp((scores + bias_for_score)*0.1)
            scores = scores*(np.exp(-1*visited_time/self.T))
            '''calculate the percentage of each score'''
            scores = scores / np.sum(scores)
            if display_state == True and node == 0:
                print('Probability distribution for current selection is %s' %scores)
            picked_node = np.random.choice(sons, size=1, p=scores)[0]
            return picked_node
        if select_mode[0] == 'KS':
            '''Kocsis and Szepesvári method from WIKI'''
            C = select_mode[1] #this is the parameter for values calculating
            args = [arg, 'visited_time_total', 'visited_time', 'activation']
            sons, res = self.GetSonAttributes(node, args)
            scores, visit_total, visit, activation = res
            score_final = scores
            '''introduce h score?'''
# =============================================================================
#             if mode_BP[0] == 'globalscore_and_ave_cx':
#                 if np.max(score_final) == 0:
#                     score_final += h_scores
# =============================================================================
            visit_final = visit_total + visit*visit_punish
            #sum_visit = max(np.sum(visit_final), 1)
            sum_visit = max(np.sum(visit_total), 1)
            if np.sum(activation) == 0: print('wrong!!!')
            '''two methods for cal values, comment one and keep one'''
            values = score_final + C * np.sqrt(
                    np.log(sum_visit) / (visit_final + 0.001)
                    )
            
            values = values * activation
            '''for root node, wo add some noise to encourage exploration'''
            '''this can be commented'''
# =============================================================================
#             if np.random.rand() > 0.3:
#                 noise = 1 - np.power(np.random.rand(len(sons)), 1)
#                 values = values * noise
# =============================================================================
            picked_index = np.argmax(values)
            picked_node = sons[picked_index]
            return picked_node, values[picked_index]
        
        if select_mode[0] == 'cx_vs_visit':
            C = select_mode[1]
            args = [arg, 'visited_time_total', 'visited_time', 'activation']
            sons, res = self.GetSonAttributes(node, args)
            scores, visit_total, visit, activation = res
            visit_final = visit_total + visit*visit_punish + 1
            score_final = (scores + 5) / visit_final
            #sum_visit = max(np.sum(visit_final), 1)
            sum_visit = max(np.sum(visit_total), 1)
            if np.sum(activation) == 0: print('wrong!!!')
            '''two methods for cal values, comment one and keep one'''
            values = score_final + C * np.sqrt(
                    np.log(sum_visit) / (visit_final + 0.001)
                    )
            
            values = (values + 10) * activation
            '''for root node, wo add some noise to encourage exploration'''
            '''this can be commented'''
# =============================================================================
#             if np.random.rand() > 0.3:
#                 noise = 1 - np.power(np.random.rand(len(sons)), 1)
#                 values = values * noise
# =============================================================================
            
            picked_node = sons[np.argmax(values)]
            return picked_node, scores[np.argmax(values)]            

        if select_mode[0] == 'KS_random':
            '''Kocsis and Szepesvári method from WIKI plan randomness'''
            power_random = select_mode[1][1]
            C = select_mode[1][0] #this is the parameter for values calculating
            args = [arg, 'visited_time_total', 'visited_time',
                         'heuristic_score', 'activation']
            sons, res = self.GetSonAttributes(node, args)
            scores, visit_total, visit, h_scores, activation = res
            score_final = scores
            visit_final = visit_total + visit*visit_punish
            #sum_visit = max(np.sum(visit_final), 1)
            sum_visit = max(np.sum(visit_total), 1)
            if np.sum(activation) == 0: print('wrong!!!')
            '''two methods for cal values, comment one and keep one'''
            values = score_final + C * np.sqrt(
                    np.log(sum_visit) / (visit_final + 0.001)
                    )
            
            values = (values + 10) * activation
            random_numbers = np.random.rand(len(sons))
            random_numbers = 1 - np.power(random_numbers, power_random)
            values = values * random_numbers
            picked_node = sons[np.argmax(values)]
            return picked_node, scores[np.argmax(values)]
        
        if select_mode[0] == 'ZHOU':
            '''use h_score'''
            bias, T = select_mode[1] #this is the parameter for values calculating
            args = [arg, 'visited_time_total', 'visited_time',
                         'heuristic_score', 'activation']
            sons, res = self.GetSonAttributes(node, args)
            scores, visit_total, visit, h_scores, activation = res
            score_final = scores + h_scores/2 + bias
            visit_final = visit_total + visit*visit_punish
            if np.sum(activation) == 0: print('wrong!!!')

            values = np.exp(-1*T*visit_final) * score_final          
            values = values * activation
            
            picked_node = sons[np.argmax(values)]
            return picked_node, scores[np.argmax(values)]
  
        if select_mode[0] == 'ZHOU2':
            '''use h_score2'''
            T_h, T_v = select_mode[1] #this is the parameter for values calculating
            args = [arg, 'visited_time_total', 'visited_time',
                         'heuristic_cost', 'activation']
            sons, res = self.GetSonAttributes(node, args)
            scores, visit_total, visit, h_costs, activation = res
            visit_final = visit_total + visit*visit_punish
            h_cost_final = np.exp(-1*T_h*visit_final) * (np.max(h_costs) - h_costs)
            score_final = scores + h_cost_final# + np.max(h_cost_final)
            score_final = score_final# + np.min(score_final) + 1
            values = np.exp(-1*T_v*visit_final) * score_final
            
            picked_node = sons[np.argmax(values)]
            return picked_node, scores[np.argmax(values)]  

        if select_mode[0] == 'ZHOU3':
            '''use h_score2'''
            bias, T_v = select_mode[1] #this is the parameter for values calculating
            args = [arg, 'visited_time_total', 'visited_time', 'activation']
            sons, res = self.GetSonAttributes(node, args)
            scores, visit_total, visit, activation = res
            visit_final = visit_total + visit*visit_punish
            score_final = scores + bias
            values = np.exp(-1*T_v*visit_final) * score_final
            
            picked_node = sons[np.argmax(values)]
            return picked_node, scores[np.argmax(values)]           
        raise(Exception('%s is not a valid selection method' %select_mode))       

    def PickBestSonNode(self, node, arg):   
        sons, (scores, visited_time)  = self.GetSonAttributes(node, 
                                                          [arg,
                                                           'visited_time_total'
                                                           ])
        if _select_mode_[0] == 'cx_vs_visit': scores = scores / visited_time
        best_son = sons[np.argmax(scores)]
        '''
        here we transcribe some data for future analysis, it won't affect results
        it can be commented to save time
        '''
        score_add = [np.min(scores), np.average(scores), np.max(scores)]
        
        '''
        if we are not so confidence, then do more selection
        These codes can be commented
        '''
# =============================================================================
#         kk = (score_add[2] / score_add[1]) - 1
#         kkk = 0
#         if kk < 0.4:
#             kkk = int(35 + (0.4 - kk) * 50)
#             self.select_mode[1] = self.select_mode[1] / 2
#         #if score_add[2] < score_add[1] * 1.4:
#             #print('unconfident decision process')
#             #print('score before', score_add[2])
#             for _ in range(kkk):
#                 self.Selection('global_score')
#             scores, visited_time, sons = self.GetSonScores(node, arg)
#             score_add = [np.min(scores), np.average(scores), np.max(scores)]
#             self.select_mode[1] = self.select_mode[1] * 2
#             #print('score after', score_add[2])
# =============================================================================
        
        '''if all son nodes' global score are 0'''
        if np.max(scores) == 0:
            print('WARNING: all scores of candidates are 0!')
            return None
        '''transcribe '''
        if log_data == True:
            self.log['score'].append(score_add)
            self.log['visit'].append([np.min(visited_time),
                                        np.average(visited_time),
                                        np.max(visited_time)])
            self.log['swap_vs_cx'].append(self.nodes[best_son]['num_total_add_gates']/\
                                          len(self.nodes[best_son]['executed_vertex']))
            self.log['num_remain_cx'].append(self.nodes[best_son]['num_remain_vertex'])
        
# =============================================================================
#         if np.argmax(scores) != np.argmax(visited_time):
#             print(np.argmax(scores), np.argmax(visited_time))
#             print(scores)
#             print(visited_time)
#             #return None
# =============================================================================
        return best_son

    def FindBestScoreInAllSons(self, node, arg):
        scores, visited_time, sons = self.GetSonScores(node, arg) 
        return np.max(scores), sons[np.argmax(scores)]
    
    def BackPropagationSingleValue(self, start_node, value, name, mode_BP):
        '''
        renew a variable reversely
        start_node: is the node the the original value extracted from, note that
                    the first backpropagated node is its father node
        value: Propagated value
        arg_name: the name (string) of this variable
        mode(string):
            '>' -> when new value > old one, renwe the old one
            '<' -> when new value < old one, renwe the old one
            'delta' -> old one = old one + (new value * args[0]^(distance))
            'globalscore' ->
                when going to a new node, we compare new_value (global score of 
                its son node)*args[0] + score (of current new node) with old 
                global score varible, if the former is larger than the latter,
                then we update the global score of this new node
            'sim_fix_CX_random' ->
        '''
        mode, args= mode_BP
        flag = True
        new_value = value
        if mode == '<':
            current_node = self.nodes[start_node]['father_node']
            while flag == True and current_node != self.root_node:
                old_value = self.nodes[current_node][name]
                if new_value < old_value:
                    self.nodes[current_node][name] = new_value
                    current_node = self.nodes[current_node]['father_node']
                    old_value = self.nodes[current_node][name]
                else:
                    flag = False
        if mode == '>':
            current_node = self.nodes[start_node]['father_node']
            while flag == True and current_node != self.root_node:
                old_value = self.nodes[current_node][name]
                if new_value > old_value:
                    self.nodes[current_node][name] = new_value
                    current_node = self.nodes[current_node]['father_node']
                else:
                    flag = False
        if mode == 'delta':
            current_node = self.nodes[start_node]['father_node']
            K = args[0]
            while current_node != self.root_node:
                old_value = self.nodes[current_node][name]
                #old_value = old_value + (new_value - old_value)*K
                old_value = old_value + new_value*K
                self.nodes[current_node][name] = old_value
                '''go to next node'''
                current_node = self.nodes[current_node]['father_node']
                K = K * K
        if mode == 'globalscore':
            current_node = self.nodes[start_node]['father_node']
            K = args[0]#recommended value is 0.7
            #print(K_list)
            while flag == True and current_node != self.root_node:
                '''calculate decay rate for propagated score'''          
                old_value = self.nodes[current_node]['global_score']
                #print(K)
                new_value = new_value * K
                new_value = self.nodes[current_node]['score'] + new_value
                if new_value > old_value:
                    self.nodes[current_node]['global_score'] = new_value
                    current_node = self.nodes[current_node]['father_node']
                else:
                    flag = False
            '''renew visited time (punishment)'''
# =============================================================================
#             while current_node != self.root_node:
#                 self.nodes[current_node]['visited_time'] += 1
#                 '''go to next node'''
#                 current_node = self.nodes[current_node]['father_node']
# =============================================================================
        if mode == 'globalscore_modified':
            current_node = self.nodes[start_node]['father_node']
            pre_node = start_node #previous node used for cal swap deacy
            #K = args[0]#recommended value is 0.7
            while flag == True and current_node != self.root_node:
                '''calculate decay rate for propagated score'''
                swap_decay = self.GetSwapDecay(pre_node)
                old_value = self.nodes[current_node]['global_score']
                new_value = new_value * swap_decay
                new_value = self.nodes[current_node]['score'] + new_value
                if new_value > old_value:
                    self.nodes[current_node]['global_score'] = new_value
                    '''go to next node'''
                    pre_node = current_node
                    current_node = self.nodes[current_node]['father_node']
                else:
                    flag = False
            '''renew visited time (punishment)'''
            while current_node != self.root_node:
                self.nodes[current_node]['visited_time'] += 1
                '''go to next node'''                
                current_node = self.nodes[current_node]['father_node']
        if mode == 'globalscore_and_ave_cx':
            amplify_num = mode_BP[1][0]
            #cal global score for start node
            num_swap = 1
            num_cx = self.nodes[start_node]['score']
            new_score = np.power(num_cx, amplify_num)/num_swap
            self.nodes[start_node]['global_score'] += new_score
            current_node = self.nodes[start_node]['father_node']
            #K = args[0]#recommended value is 0.7
            #print(K_list)
            while flag == True and current_node != self.root_node:
                '''calculate decay rate for propagated score'''          
                old_score = self.nodes[current_node]['global_score']
                num_cx += self.nodes[current_node]['score']
                num_swap += 1
                new_score = np.power(num_cx, amplify_num)/num_swap
                if new_score > old_score:
                    self.nodes[current_node]['global_score'] = new_score
                    current_node = self.nodes[current_node]['father_node']
                else:
                    flag = False
            '''renew visited time (punishment)'''
            while current_node != self.root_node:
                self.nodes[current_node]['visited_time'] += 1
                '''go to next node'''
                current_node = self.nodes[current_node]['father_node']

        if mode == 'globalscore_depth':
            K = args[0]#recommended value is 0.7
            new_value_depth = new_value
            current_node = self.nodes[start_node]['father_node']
            current_depth = self.GetLastSearchDpeth()
            #depth_threshold = np.average(self.search_depth) * args[1]
            depth_ave = self.ave_depth
            depth_min = self.min_depth
            while True and current_node != self.root_node:
                '''calculate decay rate for propagated score'''          
                old_value = self.nodes[current_node]['global_score']
                K_depth = K
                if current_depth < depth_min:
                    #print('dada')
                    K_depth = 0.99#default 0.9
                else:
                    if current_depth < depth_ave * 1:
                        K_depth = 0.95
# =============================================================================
#                     else:
#                         if current_depth < depth_ave * 0.4:
#                             K_depth = 0.8                   
# =============================================================================
                new_value = new_value * K + self.nodes[current_node]['score']
                new_value_depth = new_value_depth * K_depth + self.nodes[current_node]['score']
                '''update global value'''
                if new_value > old_value:
                    self.nodes[current_node]['global_score'] = new_value
                next_node = self.nodes[current_node]['father_node']
                '''update global_score with depth'''
                if next_node == self.root_node:
                    if new_value_depth > self.nodes[current_node]['global_score']:
                        self.nodes[current_node]['global_score'] = new_value_depth
                    #print(current_depth)
                    break
                else:
                    current_node = next_node
                    current_depth -= 1
        if mode == 'sim_fix_CX_random':
            num_swap = self.mode_sim[0]
            # store simulation result to the start_node
            self.nodes[start_node]['swap_cx_list'] = new_value
            self.nodes[start_node]['num_cx_sim'] = np.sum(new_value)
             
            new_list = new_value.copy()
            new_list.pop()
            new_list.insert(0, self.nodes[start_node]['score'])
            new_value = np.sum(new_list)
            
            current_node = self.nodes[start_node]['father_node']
            flag = True
            while flag == True and current_node != self.root_node:
                '''calculate decay rate for propagated score'''          
                old_value = self.nodes[current_node]['num_cx_sim']
                
                #num_exe_cx = self.nodes[current_node]['score']
                #new_value = new_value - num_exe_cx * new_value / num_cx + 1
                
                if new_value > old_value:
                    # store better simulation result
                    self.nodes[current_node]['swap_cx_list'] = new_list.copy()
                    self.nodes[current_node]['num_cx_sim'] = new_value
                    # update value
                    new_list.pop()
                    new_list.insert(0, self.nodes[current_node]['score'])
                    new_value = np.sum(new_list)
                    # go to next node
                    current_node = self.nodes[current_node]['father_node']
                else:
                    flag = False
# =============================================================================
#             if current_node == self.root_node:
#                 print('to root node')
#             else:
#                 print('not root node')
# =============================================================================

    def Selection(self, arg='global_score'):
        '''
        one playout, arg is decision parameter (string)
        we pick son node (randomly?) from root node to a leaf node, then, expand 
        the leaf node and back propagate the best score of expanded nodes
        return:
            added nodes and search depth
        '''
        current_node = self.root_node
        search_depth = 0
        flag_sim = True
        first_selec_node = None
        '''proceed until find leaf node'''
        while self.out_degree[current_node] != 0:
            search_depth += 1
            current_node, current_score = self.PickBestSonNodeRandomly(current_node, arg)
            # judge whether wo should ban the expanded nodes from simulation
            if current_node in self.simulated_nodes: flag_sim = False
            '''the selected node should be deactivated'''
            '''this can be commented if you want'''
# =============================================================================
#             if current_score > min_deactivation_score:
#                 self.ResetActivation(mode=['son', father_node])
#                 self.DeactivateNode(current_node)
# =============================================================================
            
            # we remember the first selected node
            if search_depth == 1:
                first_selec_node = current_node
                #first_selec_score = current_score
            '''update total visited time'''
            self.nodes[current_node]['visited_time_total'] += 1
            '''each time, we lock the node with too many visit'''
# =============================================================================
#             if self.nodes[current_node]['visited_time_total'] >= max_visit:
#                 self.DeactivateNode(current_node)
# =============================================================================
         
        '''record depth for this selection'''
        '''this can be commented'''
# =============================================================================
#         self.UpdateSearchDpeth(search_depth, first_selec_node)
# =============================================================================
        
        '''judge whether the last selected node should be deactivated'''
        '''this can be commented if you want'''
# =============================================================================
#         if search_depth >= max_depth:
#             self.DeactivateNode(current_node)
# =============================================================================
        '''judge whether the first selected node should be deactivated'''
        '''this can be commented if you want'''
# =============================================================================
#         if first_selec_score > min_deactivation_score:
#             '''each time, we lock the best son of root node to encourage exploration'''
#             self.ResetActivation(mode=['son of root'])
#             self.DeactivateNode(first_selec_node)
# =============================================================================

        '''is leaf node finished?'''
        # current_node is the expanded node
        if self.nodes[current_node]['num_remain_vertex'] != 0:
            '''If it is not finished, then expand the leaf node'''
            self.ExpandNode(current_node)
            '''back propagate the best score of newly added leaf (son) nodes'''
            best_score, best_son = self.FindBestScoreInAllSons(current_node,
                                                               arg)
            '''count heuristic score?'''
            # we may comment this because it may be moved to the h cost function
# =============================================================================
#             if _select_mode_[0] == 'ZHOU3':
#                 best_score += self.nodes[best_son]['heuristic_score']
# =============================================================================
            
            self.BackPropagationSingleValue(best_son, best_score, 'score',
                                            mode_BP)
        else:
            best_score = self.nodes[current_node][arg]
            self.BackPropagationSingleValue(current_node, best_score, 'score',
                                            mode_BP)
            #self.DeactivateNode(self.nodes[current_node]['father_node'])
# =============================================================================
#         '''if no new CNOT gate is executed, then put some punishment'''
#         self.BackPropagationSingleValue(current_node, 1, 'visited_time',
#                                         'delta', [0.8])
# =============================================================================
        return search_depth, current_node, flag_sim, first_selec_node

    def SelectionWithSimulation(self, arg):
        '''
        one playout, arg is decision parameter (string)
        after expansion, we do simultaions from best son node
        return:
            added nodes and search depth
        '''
        current_node = self.root_node
        search_depth = 0
        '''proceed until find leaf node'''
        while self.out_degree[current_node] != 0:
            search_depth += 1
            current_node = self.PickBestSonNodeRandomly(current_node, arg)
            '''update visited time'''
            self.nodes[current_node]['visited_time_total'] += 1
        '''is leaf node finished?'''
        if self.nodes[current_node]['num_remain_vertex'] != 0:
            '''Expand the leaf node'''
            self.ExpandNode(current_node)
            sim_node = current_node
            '''do simulation'''
            best_score = self.Simultation(sim_node, self.mode_sim)
            '''back propagate the best score of newly added leaf (son) nodes'''
            if best_score == None:
                best_score, _ = self.FindBestScoreInAllSons(current_node, arg)
            self.BackPropagationSingleValue(current_node, best_score, 'score',
                                            'globalscore', [score_decay_rate])
            
        return search_depth

    def DeleteNodes(self, nodes):
        '''delete nodes and all its successors'''
        '''
        WARNING: it is not finished yet because we only delete its son nodes
        right now
        '''
        #self.remove_nodes_from(nodes) # this is the old way
        for node in nodes:
            '''update son_nodes of its father node'''
            father_node = self.nodes[node]['father_node']
            if father_node != None:
                self.nodes[father_node]['son_nodes'].remove(node)
            '''delete'''
            T_succ = nx.dfs_tree(self, node)
            self.remove_nodes_from(T_succ.nodes)
            self.node_size -= len(T_succ.nodes)

    def FindCXMinDis(self,node):
        '''
        Find the executable CX in the front layer with minimum distance and
        output the CX and this distance
        '''
        executable_vertex = self.nodes[node]['executable_vertex']
        mapping = self.nodes[node]['mapping']
        min_CX_dis = 1000
        num_CX = 1
        chosen_v_DG = None
        for v in executable_vertex:
            CX = self.DG.nodes[v]['operation'].involve_qubits_list
            CX_phy = mapping.LogToPhy(CX[0]), mapping.LogToPhy(CX[1])
            currant_CX_dis = self.shortest_length_AG[CX_phy[0]][CX_phy[1]]
            if currant_CX_dis < min_CX_dis:
                min_CX_dis = currant_CX_dis
                chosen_v_DG = v
        while chosen_v_DG != None and self.DG.out_degree(chosen_v_DG) == 1:
            chosen_v_DG = list(self.DG.successors(chosen_v_DG))[0]
            num_CX += 1
        min_num_swap = min_CX_dis - 1
        if min_CX_dis == 999: min_num_swap = 0
        return min_num_swap, num_CX
    
    def Fallback(self):
        print('Fallback!')
        start_node = self.root_node
        '''find the initial node for fallback'''
        while self.nodes[start_node]['score'] == 0:
            deleted_node = start_node
            start_node = self.nodes[start_node]['father_node']
        '''extract swaps list'''
        executable_vertex = self.nodes[start_node]['executable_vertex']
        mapping = self.nodes[start_node]['mapping']
        min_CX_dis = 1000
        for v in executable_vertex:
            CX = self.DG.nodes[v]['operation'].involve_qubits_list
            CX_phy = mapping.LogToPhy(CX[0]), mapping.LogToPhy(CX[1])
            currant_CX_dis = self.shortest_length_AG[CX_phy[0]][CX_phy[1]]
            if currant_CX_dis < min_CX_dis:
                min_CX_dis = currant_CX_dis
                chosen_CX_phy = CX_phy
        path = self.shortest_path_AG[chosen_CX_phy[0]][chosen_CX_phy[1]].copy()
        #num_swap = int(np.ceil(min_CX_dis/2))
        num_swap = int(min_CX_dis - 1)
        '''set new root node and delete redunant nodes'''
        self.root_node = start_node
        self.DeleteNodes([deleted_node])
        self.nodes[start_node]['son_nodes'] = []
        '''add swaps'''
        flag = True
        for i in range(num_swap):
            if flag == True:
                added_swap = path.pop(0), path[0]
            else:
                added_swap = path.pop(), path[-1]
            flag = not flag
            added_node = self.ExpandNodeViaSWAP(self.root_node, [added_swap])
            self.root_node = added_node
        if self.nodes[self.root_node]['score'] == 0:
            '''
            if the newly added node still can't execute any CX, there is sth
            wrong with the fallback procedure
            '''
            raise(Exception('Fallback error!'))

    def FallbackOneSwap(self):
        '''UNFINISHED'''
        start_node = self.root_node
        '''extract swaps list'''
        executable_vertex = self.nodes[start_node]['executable_vertex']
        mapping = self.nodes[start_node]['mapping']
        min_CX_dis = 1000
        for v in executable_vertex:
            CX = self.DG.nodes[v]['operation'].involve_qubits_list
            CX_phy = mapping.LogToPhy(CX[0]), mapping.LogToPhy(CX[1])
            currant_CX_dis = self.shortest_length_AG[CX_phy[0]][CX_phy[1]]
            if currant_CX_dis < min_CX_dis:
                min_CX_dis = currant_CX_dis
                chosen_CX_phy = CX_phy
        path = self.shortest_path_AG[chosen_CX_phy[0]][chosen_CX_phy[1]].copy()

    
    def Decision(self):
        '''
        choose one leaf node, delete all other leaf nodes of its father node
        '''
        mode = mode_decision[0]
        father_node = self.root_node
        
        if mode == 'global_score':
            best_son = self.PickBestSonNode(father_node, 'global_score')
        
        if mode == 'sim_fix_CX':
            decay_score, decay_sim = mode_decision[1], mode_decision[2]
            sons, res  = self.GetSonAttributes(father_node, 
                                               ['global_score',
                                                'num_cx_sim'])
            score, num_cx = res
            num_swap = self.mode_sim[0]
            #num_cx = np.sum(swap_list)
            ave_cx = num_cx / (num_swap)
            final_score = decay_score * score + decay_sim * ave_cx
            best_son = sons[np.argmax(final_score)]

            '''transcribe data'''
            if log_data == True:
                self.log['num_remain_cx'].append(self.nodes[best_son]['num_remain_vertex'])
            
        if mode == 'num_sim_swap':
            #best_son = self.PickBestSonNode(father_node, 'sim_score')
            sons, (sim_score, global_score)  = self.GetSonAttributes(father_node, 
                                                                      ['sim_score',
                                                                       'global_score' ])
            #normalize
            #print('')
            #print('sim_score', sim_score)
            #print('global_score', global_score)
            #print('max sim', np.argmax(sim_score), np.max(sim_score))
            #print('max global_score', np.argmax(global_score), np.max(global_score))
            sim_score = sim_score / (np.max(sim_score)+0.001)
            global_score = global_score / (np.max(global_score)+0.001)
            chosen_index = np.argmax(sim_score + global_score)
            #print(chosen_index)
            best_son = sons[chosen_index]
            #print('best_son is', best_son)
        # if we can't find a son node, we do not upated root node
        if best_son == None: return self.root_node
        
        '''if depth is too shallow, we don't update root node'''
        '''this section can be commented'''
# =============================================================================
#         next_depth = np.average(self.nodes[best_son]['search_depth'])
#         if next_depth < self.ave_depth * 0.95 and self.first_finish_node == None:
#             print(np.average(self.nodes[best_son]['search_depth']))
#             self.nodes[best_son]['search_depth'] = []
#             return self.root_node
# =============================================================================
            
        '''reset activation'''
        self.ResetActivation()
        if delete_trivival_swap_decision == 1:
            '''
            check wehther current root_node does not exe any CX gates, if yes,
            we delete the swap that the father of root_node does to avoid dead
            loop
            '''
            if father_node != 0:
                if self.nodes[father_node]['executed_vertex_current'] == []:
                    swap_delete = self.nodes[father_node]['added_SWAP']
                    if self.nodes[best_son]['added_SWAP'] == swap_delete:
                        '''delete repeated swap'''
                        self.DeleteNodes([best_son])
                        best_son = self.PickBestSonNode(father_node,
                                                        'global_score')
        '''delete residual nodes'''
        deleted_nodes = self.nodes[father_node]['son_nodes'].copy()
        deleted_nodes.remove(best_son)
        self.DeleteNodes(deleted_nodes)
        #self.nodes[father_node]['son_nodes'] = [best_son]
        '''update fallback count'''
        if self.nodes[best_son]['score'] == 0:
            self.fallback_count += 1
        else:
            self.fallback_count = 0
        if self.fallback_count >= self.fallback_value:
            self.Fallback()
            self.fallback_count = 0
            return self.root_node
        '''update root node'''
        self.root_node = best_son
#        if display_state == True: print('Chosen next node is %d' %best_son)
        if display_state == True:
            print('\r%d gates unfinished'
                  %self.nodes[self.root_node]['num_remain_vertex'],end='')
            #print('added swap is %s\n' %self.nodes[best_son]['added_SWAP'])
# =============================================================================
#         if _select_mode_[0] == 'ZHOU2':
#             self.AddHCost(self.root_node)
# =============================================================================
        #self.PrintSonNodesArgs(self.root_node, ['heuristic_cost'])
        '''update average depth for the new root node'''
# =============================================================================
# #        print('root node ave depth before is', self.ave_depth)
#         self.ave_depth = np.average(self.nodes[self.root_node]['search_depth'])
#         self.min_depth = np.min(self.nodes[self.root_node]['search_depth'])
# #        print('root node total depth is', self.nodes[self.root_node]['search_depth'])
# #        print('root node ave depth is', self.ave_depth)
# =============================================================================
        return self.root_node
    
    def CalNodesInDGForSim(self, node, num_gates):
        executable_nodes = self.nodes[node]['executable_vertex']
        executed_nodes = self.nodes[node]['executed_vertex']
        add_nodes, _ = GetSubCxList(self.DG,
                                    executable_nodes,
                                    executed_nodes,
                                    num_gates)
        self.nodes[node]['sim_nodes'] = add_nodes
        #self.nodes[node]['num_sim_nodes'] = len(add_nodes)

    def PrintPhyCircuit(self, node=None):
        '''unfinished'''
        q_phy = QuantumRegister(20, 'v')
        cir_phy = QuantumCircuit(q_phy)
        current_node = 0
        flag = True
        CNOT_pos = 0
        gate_num = 0
        while flag == True:
            swaps = self.nodes[current_node]['added_SWAP']
            CNOT_vertex = self.nodes[current_node]['executed_vertex']
            current_mapping = self.nodes[current_node]['mapping']
            '''draw swap gates'''
            for swap in swaps:
                cir_phy.swap(q_phy[swap[0]], q_phy[swap[1]])
                #gate_num += 3
            '''draw CNOT gates'''
            if len(CNOT_vertex) > CNOT_pos:
                num_CNOT = len(CNOT_vertex) - CNOT_pos
                for i in range(num_CNOT):
                    op = self.DG.nodes[CNOT_vertex[CNOT_pos]]['operation']
                    control = current_mapping.LogToPhy(op.involve_qubits[0])
                    target = current_mapping.LogToPhy(op.involve_qubits[1])
                    cir_phy.cx(q_phy[control], q_phy[target])
                    gate_num += 1
                    CNOT_pos += 1
                cir_phy.barrier()
            '''judge whether the circuit is finished'''
            if self.nodes[current_node]['num_remain_vertex'] == 0:
                flag = False
            else:
                sons = list(self.successors(current_node))
                if len(sons) > 1: raise(Exception('more than one successors'))
                current_node = sons[0]
        fig = (cir_phy.draw(scale=0.7, filename=None, style=None, output='mpl',
                            interactive=False, line_length=None, plot_barriers=True,
                            reverse_bits=False))
        print('total CX gate number of physical circuit is %d' %gate_num)
        fig.savefig('MCTreeSearch_phy.pdf', format='pdf', papertype='a4')
        
    def PrintNodeArgs(self, node, names):
        print('  node %d' %node)
        for name in names:
            print('    %s is %s' %(name, self.nodes[node][name]))
    
    def PrintSonNodesArgs(self, father_node, names_son, names_father=[]):
        if not isinstance(names_son, list) and not isinstance(names_son, tuple):
            raise(Exception('names argument must be list or tuple, but it is %s'
                            %type(names_son)))
        if not isinstance(names_father, list) and not isinstance(names_father, tuple):
            raise(Exception('names argument must be list or tuple, but it is %s'
                            %type(names_father)))            
        print('father node is %d' %father_node)
        sons = self.nodes[father_node]['son_nodes']
        for name in names_father:
            print('    %s is %s' %(name, self.nodes[father_node][name]))
        print('all son nodes of %d' %father_node)
        for son in sons:
            self.PrintNodeArgs(son, names_son)
            
    def PrintInvolvePhyQ(self, node):
        '''print all involved physical qubits in front layer of a given node'''
        v = FindInvolvedPhyQInOneLayer(self.nodes[node]['executable_vertex'], 
                                       self.nodes[node]['mapping'], self.DG)
        print('All involved physical qubits in front layer of node %d is %s'
              %(node, v))
        
    def DrawLogData(self):
        for key in self.log.keys():
            if key == 'swap_vs_cx':
                plt.plot(self.log[key][9:])
            else:
                plt.plot(self.log[key])
            plt.legend()
            plt.grid()
            plt.title(key)
            plt.show()