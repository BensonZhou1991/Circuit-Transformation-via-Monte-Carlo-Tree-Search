# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 00:16:18 2020

@author: zxz58
"""
from circuittransform.operation import GetSubDG, ExecuteAllPossibileNodesInDG
from circuittransform.operation import SWAPInArchitectureGraph
import numpy as np

def SimultationRandomFixCX(MCTree, start_node, mode_args):
    '''
    Different modes for simulation
    ATTENTION: we won't do any Backpropagation in Simulation
    '''
    num_CX = mode_args[0]
    times_sim = mode_args[1] #simulation times
    num_sim = []
    '''get information from MCT'''
    AG = MCTree.AG
    shortest_path_AG = MCTree.shortest_path_AG
    executable_nodes = MCTree.nodes[start_node]['executable_vertex'].copy()
    executed_nodes = MCTree.nodes[start_node]['executed_vertex'].copy()
    mapping = MCTree.nodes[start_node]['mapping'].Copy()
    DG = GetSubDG(MCTree.DG, executable_nodes, executed_nodes, num_CX)
        
    num_best_swaps = 10000
    for i in range(times_sim):
        res = SimulateOneTimeFixCX(num_CX,
                                   executable_nodes.copy(),
                                   executed_nodes.copy(),
                                   mapping.Copy(), shortest_path_AG,
                                   DG, AG)
        _, num_add_swaps, swap_cx_list = res
        num_sim.append(num_add_swaps)
        if num_add_swaps < num_best_swaps:
            num_best_swaps = num_add_swaps
            swap_cx_list_best = swap_cx_list
    #print(num_sim)
    #print(np.min(num_sim))
    return swap_cx_list_best
            

def SimulateOneTimeFixCX(num_exe_CX_total,
                         executable_vertex, executed_vertex,
                         mapping, shortest_path_AG,
                         DG, AG):
    '''start a simulation from start_node executing fixed number of CNOT'''
    num_exe_CX_current = 0
    num_add_swaps = 0
    swap_cx_list = []
    #print('front layer CNOTs are', executable_vertex)

    '''execute CNOTs along shortest paths'''
    while num_exe_CX_current < num_exe_CX_total and executable_vertex != []:
        v_index = np.random.randint(len(executable_vertex))
        exe_vertex = executable_vertex[v_index]
        #print('chosen node in DG is', exe_vertex)
        #print(DG.nodes[exe_vertex]['operation'].involve_qubits_list)
        #print('current executable_vertex is', executable_vertex)
        q_c, q_t = DG.nodes[exe_vertex]['operation'].involve_qubits
        v_c, v_t = mapping.LogToPhy(q_c), mapping.LogToPhy(q_t)
        path = shortest_path_AG[v_c][v_t]
        res = ConductCNOTInDGAlongPath(DG, exe_vertex, path, mapping, AG,
                                       executable_vertex=executable_vertex,
                                       executed_vertex=executed_vertex)
        executable_vertex, executed_vertex, num_exe_CX, num_swaps, swap_cx_add = res
        #print('new executable_vertex is', executable_vertex)
        num_exe_CX_current += num_exe_CX
        num_add_swaps += num_swaps
        swap_cx_list.append(swap_cx_add)
        #print('now we have executed %d gates' %sum(CX_flag_add))
    return num_exe_CX_current, num_add_swaps, swap_cx_list

def ConductCNOTInDGAlongPath(DG, vertex, path, mapping, AG,
                             executable_vertex=None,
                             executed_vertex=None):
    '''
    conduct CNOT in a vertex in DG along a specific path([control, ..., target])
    in architecture graph RANDOMLY, then, renew physical circuit and mapping
    input:
        remove_node: will the node in DG be removed?
        exe_other_CX: should we check and execute other possible CNOT whem swapping?
                      ATTENTION: this fuction is only valid for undirected AG
    '''
    swap_cx_list = []
    num_exe_CX = 0
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
            if np.random.rand() > 0.5:
                # execute swap from another end of path next time
                flag_head = not flag_head
        else:
            SWAPInArchitectureGraph(path[v_t_pos], path[v_t_pos-1], mapping,
                                    q_phy=None, cir_phy=None, draw=False)
            v_t_pos -= 1
            if np.random.rand() > 0.5:
                # execute swap from another end of path next time
                flag_head = not flag_head
        res = ExecuteAllPossibileNodesInDG(executable_vertex, executed_vertex,
                                           AG, DG, mapping, draw=False, DiG=None,
                                           edges_DiG=None, cir_phy=None,
                                           q_phy=None, out_removed_nodes=True)
        executed_vertex, executable_vertex, removed_nodes = res
        num_exe_CX += len(removed_nodes)
        swap_cx_list.append(len(removed_nodes))
        
    return executable_vertex, executed_vertex, num_exe_CX, num_swaps, swap_cx_list