# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 00:45:34 2019

@author: zxz58

This module is for functions on operations
"""

import networkx as nx
import numpy as np

def OperationToDependencyGraph(operations):
    '''
    create dependency graph
    input:
        operations a list of all operations instances
    '''
    first_gates = [-1]*50 #a list in which each element is the node that takes the first place of qubits
    
    num_vertex = len(operations)
    DG = nx.DiGraph()
    num_q_log = 0
    DG.add_nodes_from(list(range(num_vertex)))
    for i in range(num_vertex):
        current_operation = operations[i]
        qubits = current_operation.InvolveQubitsList()
        for qubit in qubits:
            if first_gates[qubit] == -1: first_gates[qubit] = i
            if qubit + 1 > num_q_log: num_q_log = qubit + 1
        DG.add_node(i, operation = current_operation)
        if current_operation.dependent_operations != []:
            DG.add_node(i, root = False)
            for current_de in current_operation.dependent_operations:
                DG.add_edge(operations.index(current_de), i)
        else:
            DG.add_node(i, root = True)
    DG.first_gates = first_gates
    DG.num_q_log = num_q_log
    return DG

def FindExecutableNode(dependency_graph,
                       executed_vertex=None,
                       executable_vertex=None,
                       removed_vertexes=None):
    '''
    WHEN executable_vertex = None:
        Use dependency graph to find the executable vertexes/nodes, i.e., nodes
        in current level
        return:
            executable_nodes: a list of nodes. If no executable node, return []
    WHEN both executed_vertex and executable_vertex != None:
        only update executed_vertex and executable_vertex according to newly
        executed gates (removed_vertexes)
        return:
            executable_vertex      
    '''
    DG = dependency_graph
    if executable_vertex == None:
        degree = DG.in_degree
        executable_nodes = []
        for i in degree:
            if i[1] == 0:
                executable_nodes.append(i[0])
    else:
        executable_nodes = executable_vertex
        for removed_vertex in removed_vertexes:
            if not removed_vertex in executable_vertex: raise Exception('removed node is not executable')
            candidate_nodes = DG.successors(removed_vertex)
            executable_nodes.remove(removed_vertex)
            executed_vertex.append(removed_vertex)
            #DG.remove_node(removed_vertex)
            for node in candidate_nodes:
                flag_add = True
                '''check whether this node is executable'''
                for pre_node in DG.predecessors(node):
                    if not pre_node in executed_vertex:
                        flag_add = False
                        break
                if flag_add == True: executable_nodes.append(node)
    return executable_nodes
    
def FindExecutableOperation(dependency_graph, executable_node = None):
    '''
    Use dependency graph to find the executable operations
    return:
        executable_operation: a list of operations. If no executable operation, return []
    '''
    DG = dependency_graph
    executable_operation = []
    if executable_node == None:
        executable_node = FindExecutableNode(DG)
    for node in executable_node:
        executable_operation.append(DG.node[node]['operation'])
    return executable_operation

def ConductOperationInVertex(DG, vertex, mapping, cir_phy, q_phy, remove_node=True):
    '''
    Conduct operation in physical quantum circuit represented by
    a node/vertex in dependency graph(directed graph), then, erase the corresponding
    node/vertex in dependency graph
    '''
    conduct_operation = DG.node[vertex]['operation']
    q_c = conduct_operation.control_qubit
    q_t = conduct_operation.target_qubit
    v_c = mapping.DomToCod(q_c)
    v_t = mapping.DomToCod(q_t)
    q_phy_c = q_phy[v_c]
    q_phy_t = q_phy[v_t]
    conduct_operation.ConductOperationOutside(cir_phy, q_phy_c, q_phy_t)
    if remove_node == True: DG.remove_node(vertex)

def ConductCNOTOperationInVertex(DG, vertex, mapping, cir_phy, q_phy, reverse_drection=False, remove_node=True):
    '''
    conduct CNOT operation in physical quantum circuit represented by
    a node/vertex in dependency graph(directed graph), then, erase the corresponding
    node/vertex in dependency graph
    input:
        convert_drection -> whether use 4 H gates to reverse direction of CNOT
        remove_node -> whether the node in DG shuold be removed
    '''
    conduct_operation = DG.node[vertex]['operation']
    q_c = conduct_operation.control_qubit
    q_t = conduct_operation.target_qubit
    v_c = mapping.DomToCod(q_c)
    v_t = mapping.DomToCod(q_t)
    q_phy_c = q_phy[v_c]
    q_phy_t = q_phy[v_t]
    if reverse_drection == True:
        cir_phy.h(q_phy_c)
        cir_phy.h(q_phy_t)
        cir_phy.cx(q_phy_t, q_phy_c)
        cir_phy.h(q_phy_c)
        cir_phy.h(q_phy_t)
    else:
        cir_phy.cx(q_phy_c, q_phy_t)    
    if remove_node == True: DG.remove_node(vertex)
    
def SWAPInArchitectureGraph(vertex0, vertex1, mapping, q_phy=None, cir_phy=None, draw=True):
    '''
    Conduct SWAP in physical qubits represented by nodes in architerture graph
    Then, renew physical circuit and mapping
    '''
    if draw == True: cir_phy.swap(q_phy[vertex0], q_phy[vertex1])
    mapping.RenewMapViaExchangeCod(vertex0, vertex1)

def ConductCNOTInDGAlongPath(DG, vertex, path, mapping, draw, remove_node=True,
                             q_phy=None, cir_phy=None, edges_DiG=None,
                             exe_other_CX = False, executable_vertex=None,
                             executed_vertex=None, G=None):
    '''
    conduct CNOT in a vertex in DG along a specific path([control, ..., target])
    in architecture graph, then, renew physical circuit and mapping
    input:
        remove_node: will the node in DG be removed?
        exe_other_CX: should we check and execute other possible CNOT whem swapping?
                      ATTENTION: this fuction is only valid for undirected AG
    '''
    add_gates_count = 0
    v_c_pos = 0
    v_t_pos = len(path) - 1
    num_swaps = len(path) - 2
    flag_head = True #decide in which way (control to target or target to control) we perform SWAP
    if edges_DiG == None:
        for i in range(num_swaps):
            add_gates_count += 3
            if flag_head == True:
                SWAPInArchitectureGraph(path[v_c_pos], path[v_c_pos+1], mapping,
                                        q_phy, cir_phy, draw)
                v_c_pos += 1
                flag_head = not flag_head
            else:
                SWAPInArchitectureGraph(path[v_t_pos], path[v_t_pos-1], mapping,
                                        q_phy, cir_phy, draw)
                v_t_pos -= 1
                flag_head = not flag_head
            if exe_other_CX == True:
                res = ExecuteAllPossibileNodesInDG(executable_vertex, executed_vertex,
                                                   G, DG, mapping, draw=False, DiG=None,
                                                   edges_DiG=None, cir_phy=None,
                                                   q_phy=None)
                executable_vertex, executed_vertex = res
        if draw == True:
            ConductCNOTOperationInVertex(DG, vertex, mapping, cir_phy, q_phy,
                                         reverse_drection=False,
                                         remove_node=remove_node)
            cir_phy.barrier()
        else:
            if remove_node == True: DG.remove_node(vertex)
    else:
        for i in range(num_swaps):
            add_gates_count += 7
            if not ((path[v_c_pos], path[v_c_pos+1]) in edges_DiG) == True:
                SWAPInArchitectureGraph(path[v_c_pos], path[v_c_pos+1], mapping,
                                        q_phy, cir_phy, draw)
                v_c_pos += 1
            else:
                SWAPInArchitectureGraph(path[v_t_pos], path[v_t_pos-1], mapping,
                                        q_phy, cir_phy, draw)
                v_t_pos -= 1
        flag_4H = not ((path[v_c_pos], path[v_t_pos]) in edges_DiG)
        add_gates_count += flag_4H * 4
        if draw == True:
            ConductCNOTOperationInVertex(DG, vertex, mapping, cir_phy, q_phy,
                                         flag_4H, remove_node)
            cir_phy.barrier()
        else:
            if remove_node == True: DG.remove_node(vertex)   
    if exe_other_CX == False:
        return add_gates_count
    else:
        return add_gates_count, executable_vertex, executed_vertex
    
def RemoveUnparallelEdge(remain_edge, remove_edge):
    '''
    this function if for function FindAllPossibleSWAPParallel
    given a remain_edge set and a to be removed edge, this function will remove
    the given edge and others edges that are effected by this edge
    '''
    remain_edge.remove(remove_edge)
    n1 = remove_edge[0]
    n2 = remove_edge[1]
    iterate_edge = remain_edge.copy()
    for current_edge in iterate_edge:
        c1 = current_edge[0]
        c2 = current_edge[1]
        if (n1 == c1) or (n1 == c2) or (n2 == c1) or (n2 == c2):
            remain_edge.remove(current_edge)

def AddEdgeToList(total_swap, remain_edge):
    '''
    this is only for function FindAllPossibleSWAPParallel
    '''
    basic_combination = total_swap[-1].copy()
    for new_edge in remain_edge:
        #print(remain_edge)  
        new_combination = basic_combination.copy()
        new_combination.append(new_edge)
        total_swap.append(new_combination)
        next_remain_edge = remain_edge.copy()
        RemoveUnparallelEdge(next_remain_edge, new_edge)
        #print(next_remain_edge)
        if next_remain_edge != []:
            AddEdgeToList(total_swap, next_remain_edge)
        #print(remain_edge)

def RemoveRepetitiveSWAPCombination(total_swap):
    '''
    this is only for function FindAllPossibleSWAPParallel
    '''
    
    '''sort'''
    sort_key = lambda arg: arg[0]
    for current_swap in total_swap:
        current_swap.sort(key=sort_key)
        
    '''delete repetitive SWAP combinations'''
    i = -1
    while i < (len(total_swap) - 2):
        i += 1
        j = i
        while j < (len(total_swap) - 1):
            j += 1
            if total_swap[i] == total_swap[j]:
                total_swap.pop(j)
                j -= 1
    
def FindAllPossibleSWAPParallel(G, availiavle_vertex=None):
    '''
    find all possible combinations of SWAPs that can conducted in parallel
    return:
        list of vertex pairs each representing a SWAP, i.e., (((v00, v01), (v10, v11)...)...)
    '''
    total_swap = []
    if availiavle_vertex == None:
        node = list(G.node)
    else:
        node = availiavle_vertex
    edge = list(G.edges).copy()
    for current_edge in edge:
        total_swap.append([current_edge])
        remain_edge = edge.copy()
        RemoveUnparallelEdge(remain_edge, current_edge)
        if remain_edge != []:
            AddEdgeToList(total_swap, remain_edge)
    RemoveRepetitiveSWAPCombination(total_swap)
        
    return total_swap

def CalRemoteCNOTCostinArchitectureGraph(path, DiG=None, shortest_length_G=None):
    '''Calculate the number of CNOT in remote CNOT implementation, including the target CNOT operation'''
    '''if the architecture graph is directed, the result will be added the possible 4 H gates'''
    dis = len(path) - 1
    if DiG != None: edges = list(DiG.edges())
    if dis == 2:
        if DiG == None:
            CNOT_cost = 4
        else:
            CNOT_cost = 4 + CheckCNOTNeedConvertDirection(path[0], path[1], path[0:2], edges)*2*4 + \
            CheckCNOTNeedConvertDirection(path[1], path[2], path[1:3], edges)*2*4
    else:
        if dis ==3:
            if DiG == None:
                CNOT_cost = 6
            else:
                CNOT_cost = 6 + CheckCNOTNeedConvertDirection(path[0], path[1], path[0:2], edges)*2*4 + \
                CheckCNOTNeedConvertDirection(path[1], path[2], path[1:3], edges)*2*4 + \
                CheckCNOTNeedConvertDirection(path[2], path[3], path[2:4], edges)*2*4
    return CNOT_cost

def RemoteCNOTinArchitectureGraph(path, cir_phy, q_phy, DiG=None):
    '''
    implement remote CNOT in physical circuit via path of nodes in architecture graph, i.e., [v_c, ..., v_t]
    '''
    num_CNOT = 0
    if DiG != None: edges = list(DiG.edges())
    dis = len(path) - 1
    q = q_phy
    v_c = path[0]
    v_t = path[-1]
    if dis == 2:
        flag_4H_1 = False
        flag_4H_2 = False
        if DiG != None:
            if CheckCNOTNeedConvertDirection(path[0], path[1], path[0:2], edges) == True:
                flag_4H_1 = True
                num_CNOT += 2*4
            if CheckCNOTNeedConvertDirection(path[1], path[2], path[1:3], edges) == True:
                flag_4H_2 = True
                num_CNOT += 2*4
        if flag_4H_1 == True:
            cir_phy.h(q[path[0]])
            cir_phy.h(q[path[1]])
        cir_phy.cx(q[v_c], q[path[1]])
        if flag_4H_1 == True:
            cir_phy.h(q[path[0]])
            cir_phy.h(q[path[1]])
        if flag_4H_2 == True:
            cir_phy.h(q[path[1]])
            cir_phy.h(q[path[2]])            
        cir_phy.cx(q[path[1]], q[v_t])
        if flag_4H_2 == True:
            cir_phy.h(q[path[1]])
            cir_phy.h(q[path[2]])  
        if flag_4H_1 == True:
            cir_phy.h(q[path[0]])
            cir_phy.h(q[path[1]])        
        cir_phy.cx(q[v_c], q[path[1]])
        if flag_4H_1 == True:
            cir_phy.h(q[path[0]])
            cir_phy.h(q[path[1]])
        if flag_4H_2 == True:
            cir_phy.h(q[path[1]])
            cir_phy.h(q[path[2]])  
        cir_phy.cx(q[path[1]], q[v_t])
        if flag_4H_2 == True:
            cir_phy.h(q[path[1]])
            cir_phy.h(q[path[2]])  
            
        num_CNOT = num_CNOT + 4
    else:
        if dis == 3:
            flag_4H_1 = False
            flag_4H_2 = False
            flag_4H_3 = False
            if DiG != None:
                if CheckCNOTNeedConvertDirection(path[0], path[1], path[0:2], edges) == True:
                    flag_4H_1 = True
                    num_CNOT += 2*4
                if CheckCNOTNeedConvertDirection(path[1], path[2], path[1:3], edges) == True:
                    flag_4H_2 = True
                    num_CNOT += 2*4    
                if CheckCNOTNeedConvertDirection(path[2], path[3], path[2:4], edges) == True:
                    flag_4H_3 = True
                    num_CNOT += 2*4   
            if flag_4H_1 == True:
                cir_phy.h(q[path[0]])
                cir_phy.h(q[path[1]])
            cir_phy.cx(q[v_c], q[path[1]])
            if flag_4H_1 == True:
                cir_phy.h(q[path[0]])
                cir_phy.h(q[path[1]])
            if flag_4H_3 == True:
                cir_phy.h(q[path[2]])
                cir_phy.h(q[path[3]])               
            cir_phy.cx(q[path[2]], q[v_t])
            if flag_4H_3 == True:
                cir_phy.h(q[path[2]])
                cir_phy.h(q[path[3]])   
            if flag_4H_2 == True:
                cir_phy.h(q[path[1]])
                cir_phy.h(q[path[2]])            
            cir_phy.cx(q[path[1]], q[path[2]])
            if flag_4H_2 == True:
                cir_phy.h(q[path[1]])
                cir_phy.h(q[path[2]])   
            if flag_4H_1 == True:
                cir_phy.h(q[path[0]])
                cir_phy.h(q[path[1]])            
            cir_phy.cx(q[v_c], q[path[1]]) 
            if flag_4H_1 == True:
                cir_phy.h(q[path[0]])
                cir_phy.h(q[path[1]])
            if flag_4H_2 == True:
                cir_phy.h(q[path[1]])
                cir_phy.h(q[path[2]])               
            cir_phy.cx(q[path[1]], q[path[2]])
            if flag_4H_2 == True:
                cir_phy.h(q[path[1]])
                cir_phy.h(q[path[2]])   
            if flag_4H_3 == True:
                cir_phy.h(q[path[2]])
                cir_phy.h(q[path[3]])               
            cir_phy.cx(q[path[2]], q[v_t])
            if flag_4H_3 == True:
                cir_phy.h(q[path[2]])
                cir_phy.h(q[path[3]])               
            
            num_CNOT = num_CNOT + 6
                   
    return num_CNOT

def IsVertexInDGOperatiable(vertex, DG, G, mapping):
    '''check whether the vertex of DG can be executed, ignoring the dependency'''
    op = DG.node[vertex]['operation']
    if len(op.involve_qubits) == 1: return True #single-qubit gate
    q0 = op.involve_qubits[0]
    q1 = op.involve_qubits[1]
    v0 = mapping.DomToCod(q0)
    v1 = mapping.DomToCod(q1)
    if (v0, v1) in G.edges():
        return True
    else:
        return False
    
def CheckCNOTNeedConvertDirection(v_c, v_t, path, edges):
    l = len(path)
    for i in range(l-1):
        c_i = path[i]
        t_i = path[i+1]
        if (c_i, t_i) in edges:
            return False
    
    return True

def CheckCNOTNeedConvertDirection2(vertex, DG, mapping, edges):
    '''input is a node in DG'''
    op = DG.node[vertex]['operation']
    q0 = op.involve_qubits[0]
    q1 = op.involve_qubits[1]
    v0 = mapping.DomToCod(q0)
    v1 = mapping.DomToCod(q1)
    if (v0, v1) in edges:
        return False
    if (v1, v0) in edges:
        return True
    raise Exception('this CNOT can not be executed')

def CheckSWAPInvolved(swaps, executable_vertex, DG, mapping):
    '''
    check whether the swaps has swap having no effect on any executable gate,
    if yes, return False
    '''
    q_phy = []
    for vertex in executable_vertex:
        op = DG.node[vertex]['operation']
        q0 = op.involve_qubits[0]
        q1 = op.involve_qubits[1]
        v0 = mapping.DomToCod(q0)
        v1 = mapping.DomToCod(q1)
        q_phy.extend([v0, v1])
    for swap in swaps:
        if (not (swap[0] in q_phy)) and (not (swap[1] in q_phy)):
            return False
    return True

def ExecuteAllPossibileNodesInDG(executable_vertex, executed_vertex, G, DG,
                                 mapping, draw, DiG, edges_DiG, cir_phy=None,
                                 q_phy=None, out_removed_nodes=False):
    '''check whether this window already has appliable vertexes, if has, then execute them'''
    temp = True
    removed_nodes_all = []
    while temp == True:
        temp = False
        removed_nodes = []
        for vertex in executable_vertex :
            if IsVertexInDGOperatiable(vertex, DG, G, mapping) == True:
                '''check whether this CNOT needs 4 H gates to convert direction'''
                if DiG != None:
                    flag_4H = ct.CheckCNOTNeedConvertDirection2(vertex, DG, mapping, edges_DiG)
                    if flag_4H == False:
                        '''if no need 4 extra H, then execute it'''
                        if draw == True:
                            ct.ConductCNOTOperationInVertex(DG, vertex, mapping, cir_phy, q_phy, reverse_drection=flag_4H, remove_node=False)
                            cir_phy.barrier()
                        removed_nodes.append(vertex)
                        temp = True
                else:
                    '''if architecture graph is undirected, execute it'''
                    if draw == True:
                        ct.ConductCNOTOperationInVertex(DG, vertex, mapping, cir_phy, q_phy, reverse_drection=False, remove_node=False)
                        cir_phy.barrier()
                    removed_nodes.append(vertex)                           
                    temp = True
        if temp == True:
            removed_nodes_all.extend(removed_nodes)
            executable_vertex = FindExecutableNode(DG, executed_vertex,
                                                   executable_vertex,
                                                   removed_nodes)
    if out_removed_nodes == False:
        return executed_vertex, executable_vertex
    else:
        return executed_vertex, executable_vertex, removed_nodes_all

def FindInvolvedPhyQInOneLayer(exectable_vertex, mapping, DG):
    '''
    Given all vertexe indexes of DG in one layer and a mapping, output all invilved 
    physical qubits indexes
    '''
    v = []
    for vertex in exectable_vertex:
        qubits = DG.node[vertex]['operation'].involve_qubits
        for q in qubits:
            v.append(mapping.LogToPhy(q))
    return v

def GetNextLayerNodeInDG(executable_nodes, DG, num_q_log):
    '''
    Given front layer nodes in DG and DG,
    Return executable nodes in the next layer
    WARNING: UNFINISHED
    '''
    print('This is an unfinished function!')
    occupied_qubits = []
    qubits_next_belongs = [-1] * num_q_log #qubit by index belong to which node
    candidate_nodes = []
    for node in executable_nodes:
        '''update currant occupied logical qubits'''
        op = DG.node[node]['operation']
        involve_qubits = op.involve_qubits_list
        occupied_qubits.extend(involve_qubits)
        '''update candidate nodes for next layer'''
        for node_next in DG.successors(node):
            if not node_next in candidate_nodes:
                candidate_nodes.append(node_next)
        
def GetSubDG(DG, executable_nodes, executed_nodes, num_total_nodes):
    '''
    Get a new DG which is the sub graph of input DG. It's node number is
    no less than num_total_nodes or
    less than num_total_nodes and all nodes are removed, i.e., the whole
    circuit is finished
    and the nodes are input executable_nodes and their successors
    We won't do any change to lists of executable_nodes and executed_nodes
    '''
    num_current_nodes = 0
    add_nodes = []
    executable_nodes = executable_nodes.copy()
    executed_nodes = executed_nodes.copy()
    while num_current_nodes < num_total_nodes and executable_nodes != []:
        num_current_nodes += len(executable_nodes)
        add_nodes.extend(executable_nodes)
        removed_nodes = executable_nodes.copy()
        executable_nodes = FindExecutableNode(dependency_graph=DG,
                                              executed_vertex=executed_nodes,
                                              executable_vertex=executable_nodes,
                                              removed_vertexes=removed_nodes)
        
    sub_DG = DG.subgraph(add_nodes).copy()
    return sub_DG

def GetSubCxList(DG, executable_nodes, executed_nodes, min_num_gates):
    num_current_nodes = 0
    add_nodes = []
    executable_nodes = executable_nodes.copy()
    executed_nodes = executed_nodes.copy()
    while num_current_nodes < min_num_gates and executable_nodes != []:
        num_current_nodes += len(executable_nodes)
        add_nodes.extend(executable_nodes)
        removed_nodes = executable_nodes.copy()
        executable_nodes = FindExecutableNode(dependency_graph=DG,
                                              executed_vertex=executed_nodes,
                                              executable_vertex=executable_nodes,
                                              removed_vertexes=removed_nodes)
    return add_nodes, num_current_nodes
    