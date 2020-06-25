# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 04:16:25 2019

@author: zxz58
"""
import circuittransform as ct
from circuittransform.inputs.inputgenerator import CreateDGfromQASMfile
from circuittransform.operation import FindExecutableNode, ExecuteAllPossibileNodesInDG
from circuittransform import Map
from circuittransform.inputs.operationU import OperationSWAP
import networkx as nx

def FindLongestLength(DG, nodes_in, nodes_out):
    currant_length = 0
    currant_nodes = nodes_in.copy()
    while currant_nodes != []:
        currant_length += 1
        next_nodes = []
        for node in currant_nodes:
            nodes_suc = DG.successors(node)
            for node_suc in nodes_suc:
                if not node_suc in nodes_out and not node_suc in next_nodes:
                    next_nodes.append(node_suc)
        currant_nodes = next_nodes.copy()
    return currant_length

def CalDepthDG(DG):
    '''cal depth for a given DG'''
    nodes_in = []
    nodes_out = []
    for node in DG.nodes():
        if DG.in_degree(node) == 0:
            nodes_in.append(node)
        if DG.out_degree(node) == 0:
            nodes_out.append(node)
    depth = FindLongestLength(DG, nodes_in, nodes_out) + 1
    return depth

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

def RenewExeLists(executable_vertex, executed_vertex, AG, DG, mapping):
    res = ExecuteAllPossibileNodesInDG(executable_vertex, executed_vertex, AG, 
                                       DG, mapping,
                                       out_removed_nodes=True,
                                       draw=False, DiG=None, edges_DiG=None)
    
    return res

def AddSwapToDG(swap, DG, mapping, q_log_to_node):
    index_added_node = DG.number_of_nodes()
    edges_delete = []
    edges_add = []
    for q_phy in swap:
        '''we want to insert swap between node_1 and 2 in DG'''
        node_1 = None
        node_2 = None
        '''find node_1 and node_2'''
        #print(mapping.MapToList())
        #print(q_phy)
        q_log = mapping.PhyToLog(q_phy)
        #print(q_log)
        #q_log = q_log.index
        node = q_log_to_node[q_log]
        '''find deleted edge'''
        if node != -1:
            node_1 = node
            node_2_candidates = DG.successors(node)
            for node_suc in node_2_candidates:
                op = DG.node[node_suc]['operation']
                q_sucs = op.InvolveQubitsList()
                for q_suc in q_sucs:
                    if q_suc == q_log: node_2 = node_suc
            if node_2 != None:
                edges_delete.append((node_1, node_2))
                edges_add.append((node_1, index_added_node))
                edges_add.append((index_added_node, node_2))
            else:
                edges_add.append((node_1, index_added_node))
        else:
            for first_node in DG.first_gates:
                if first_node == -1: continue
                op = DG.node[first_node]['operation']
                q_firsts = op.InvolveQubitsList()
                for q_first in q_firsts:
                    if q_first == q_log:
                        edges_add.append((index_added_node, first_node))
    '''delete edges'''
    swap_log = mapping.PhyToLog(swap[0]), mapping.PhyToLog(swap[1])
    #print(swap_log)
    for edge in edges_delete:
        if edge in DG.edges(): DG.remove_edge(edge[0], edge[1])
    '''add node'''
    DG.add_node(index_added_node)
    op_add = OperationSWAP(swap_log[0], swap_log[1])
    DG.node[index_added_node]['operation'] = op_add
    '''add edges'''
    for edge in edges_add:
        DG.add_edge(edge[0], edge[1])
    '''update q_log_to_node'''
    q_log_to_node[swap_log[0]] = index_added_node
    q_log_to_node[swap_log[1]] = index_added_node
        
# =============================================================================
#         if node != -1:
#             if node in executable_vertex:
#                 node_2 = node
#                 node_1_candidates = DG.predecessors(node)
#                 for node_pre in node_1_candidates:
#                     op = DG.node[node_pre]['operation']
#                     q_pres = op.inovlve_qubits_list()
#                     for q_pre in q_pres:
#                         if q_pre == q_log: node_1 = node_pre
#             else:
#                 node_1 = node
#                 node_2_candidates = DG.successors(node)
#                 for node_suc in node_2_candidates:
#                     op = DG.node[node_suc]['operation']
#                     q_sucs = op.inovlve_qubits_list()
#                     for q_suc in q_sucs:
#                         if q_suc == q_log: node_2 = node_suc
#             if node_1 != None and node_2 != None:
#                 edges_delete.append((node_1, node_2))
# =============================================================================


def CalDepth(swap_list, initial_map_list, cir_name, AG):
    '''
    AG must be undirected and you should generate it using your program rather
    than mine because the index of nodes in AG maybe different between different
    program implementation
    swap_list is [(s11, s12), (s21, s22),...] where 's' is the phy qubit
    ini_map_list: is a list representing ini mapping from log to phy qubits,
                    more specifically,
                    the index is the log q and value the phy
    cir_name is the name of the circuir, e.g., 'xor5_254.qasm'
    '''
    res = CreateDGfromQASMfile(cir_name, flag_single=True)
    #cir = res[0]
    DG, num_unidentified_gates, q_log, operations = res[1]
    DG_swap = DG.copy()
    DG_swap.first_gates = DG.first_gates
    swaps = swap_list.copy()
#    print('swaps are', swaps)
    mapping = Map(list(range(AG.number_of_nodes())), AG, initial_map_list)
    
    q_log_to_node = [-1] * 20
    
    executable_vertex = FindExecutableNode(DG)
    executed_vertex = []
    res = RenewExeLists(executable_vertex, executed_vertex, AG, DG, mapping)
    executed_vertex, executable_vertex, executed_vertex_current = res
    
    while executable_vertex != []:
        '''update q_log_to_node'''
        for node in executed_vertex_current:
            op = DG.node[node]['operation']
            q_in = op.InvolveQubitsList()
            #print(q_in)
            for q in q_in:
                q_log_to_node[q] = node
        '''exe and add swap till we can exe other gates'''
        if swaps == []:
            raise(Exception('swaps are wrong!'))
        swap = swaps.pop(0)
        AddSwapToDG(swap, DG_swap, mapping, q_log_to_node)
        mapping.RenewMapViaExchangeCod(swap[0], swap[1])
        #print(mapping.MapToList())
        '''update exe lists'''
        res = RenewExeLists(executable_vertex, executed_vertex, AG, DG, mapping)
        executed_vertex, executable_vertex, executed_vertex_current = res
    
    depth_before = CalDepthDG(DG)
    depth_after = CalDepthDG(DG_swap)
    #nx.draw_networkx(DG_swap)
    #print(depth2, depth)
    return depth_before, depth_after

if __name__ == '__main__':
    AG = ct.GenerateArchitectureGraph(20, ['IBM QX20'])
    swaps = [(3,2),(1,2),(1,2),(3,2),(1,2),(4,3),(3,2),(1,2)]
    depths = CalDepth(swaps, list(range(16)), 'xor5_254.qasm', AG)
    print('depth before is %d and after %d' %depths)