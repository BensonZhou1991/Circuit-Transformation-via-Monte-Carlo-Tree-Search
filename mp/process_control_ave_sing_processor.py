# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 13:18:29 2020

@author: zxz58
"""
import sys
sys.path.append('..')
import multiprocessing as mp
import numpy as np
from to_python import to_python
import time
sim_cpp = to_python.SimTest
decay_BP = 0.7 # recommend: 0.7

def GetCxLists(DG, sim_nodes, removed_nodes):
    gate0 = []
    gate1 = []
    #print('all node', sim_nodes)
    #print('removed', removed_nodes)
    for node in sim_nodes:
        if not node in removed_nodes:
            q0, q1 = DG.nodes[node]['operation'].involve_qubits_list
            gate0.append(q0)
            gate1.append(q1)
    #print('gate0', gate0)
    return gate0, gate1

def IniMultiProcessForSim(process_num):
    return None, None, None 
        
def StartNewSim(queue_in,
                start_node,
                first_selec_node,# node that add the sim_score
                depth,
                MCTree,
                sim_num):
    DG = MCTree.DG
    #root_node = MCTree.root_node
    
    remove_nodes = []
    #update removed nodes
# =============================================================================
#     remove_nodes = []
#     current_node = start_node
#     while current_node != MCTree.root_node:
#         remove_nodes.extend(MCTree.nodes[current_node]['executed_vertex_current'])
#         current_node = MCTree.nodes[current_node]['father_node']
# =============================================================================
    
    sim_nodes = MCTree.nodes[start_node]['sim_nodes']
    gate0, gate1 = GetCxLists(DG, sim_nodes, remove_nodes)
    if len(gate0) <= 1: return False
    # gen map list for simulation
    mapping = MCTree.nodes[start_node]['mapping'].MapToList()
    for i in range(20):
        # complete map list
        if not i in mapping: mapping.append(i)
    
    num_swap_sim = sim_cpp(gate0, gate1, mapping, sim_num)
    num_gates = len(gate0)
    BackPropagatSimRes(MCTree, start_node, num_swap_sim, num_gates)
    return True
    
def GetSimRes(MCTree, queue_out,
              finish_all=False,
              sim_count = 1):
    return 0

def HaltAllSim(p_list):
    pass

def BackPropagatSimRes(MCTree, start_node, num_swap_sim, num_gates):
    if start_node == MCTree.root_node:
        return None
    
    # first way for sim_score
    #h_score = num_gates / num_swap_sim
    #h_score = h_score * decay_BP
    # second way for sim_score
    #print(num_gates, num_swap_sim)
    h_score = num_gates * np.power(decay_BP, num_swap_sim/2)
    
    #print(start_node, num_gates, num_swap_sim, h_score)
    
    # we backpropagation heuristic score
    new_value = MCTree.nodes[start_node]['score'] + h_score
    if new_value > MCTree.nodes[start_node]['global_score']:
        MCTree.nodes[start_node]['global_score'] = new_value
        MCTree.BackPropagationSingleValue(start_node,
                                          value=new_value,
                                          name=None,
                                          mode_BP=['globalscore',[decay_BP]])
    