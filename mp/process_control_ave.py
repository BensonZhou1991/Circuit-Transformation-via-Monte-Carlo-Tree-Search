#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 01:46:00 2020

@author: xianzhou
"""
import multiprocessing as mp
import numpy as np
from to_python import to_python
import time
sim_c = to_python.SimTest

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
    queue_in = mp.Queue()
    queue_out = mp.Queue()
    process_list = []
    for i in range(process_num):
        #print('start to generate the %dth processor'%i)
        p = mp.Process(target=IniProcessForSim,args=(queue_in, queue_out,))
        p.start()
        #print('finish to generate the %dth processor'%i)
        process_list.append(p)
    return process_list, queue_in, queue_out

def IniProcessForSim(queue_in, queue_out):
    while True:
        try:
            MCT_info = queue_in.get(block=False)
        except:
            pass
        else:
            #t_s = time.time()
            start_node = MCT_info[0]
            map_list = MCT_info[1]
            sim_times = MCT_info[2]
            gate0 = MCT_info[3]
            gate1 = MCT_info[4]
            
            num_gates = len(gate0)
            #print('start sim', first_selec_node)
            #print('gate', len(gate0))
            #print(gate0)
            #print(gate1)
            #print('map', map_list)
            #t_s = time.time()
            num_swap_sim = sim_c(gate0, gate1, map_list, sim_times)
            #num_swap_sim += depth
            #print('finish sim, score', num_swap_sim)
            queue_out.put((start_node, num_swap_sim, num_gates))
            #print('time for sim is ', time.time()-t_s)
        
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
    
    mapping = MCTree.nodes[start_node]['mapping'].MapToList()
    for i in range(20):
        # complete map list
        if not i in mapping: mapping.append(i)
    
    #print('sim put')
    #print('gate', len(gate0))
    #print('map', mapping)
    queue_in.put((start_node,
                  mapping, sim_num,
                  gate0, gate1))
    return True
    
def GetSimRes(MCTree, queue_out,
              finish_all=False,
              sim_count = 1):
    '''obtain sim result and BP the score'''
    #if finish_all == True and sim_count != 0:
    #    print('selec finished but still %d sim remains' %sim_count)
    finish_sim_count = 0
    #if finish_all == True: print('wait')
    #if finish_all == True: print('sim_count', sim_count)
    while sim_count >= 1:
        try:
            start_node, num_swap_sim, num_gates = queue_out.get(block=finish_all)
        except:
            break
        else:
            finish_sim_count += 1
            #num_gates = MCTree.nodes[MCTree.root_node]['num_sim_nodes']
            BackPropagatSimRes(MCTree, start_node, num_swap_sim, num_gates)
            if finish_all == True:
                sim_count -= 1
                #if finish_all == True: print('sim_count', sim_count)
    return finish_sim_count


def HaltAllSim(p_list):
    for p in p_list:
        #p.join()
        p.terminate()
        
# =============================================================================
# def BackPropagatSimRes(MCTree, first_selec_node, num_swap_sim, num_gates):
#     # we backpropagation heuristic score
#     current_score = MCTree.nodes[first_selec_node]['sim_score']
#     new_score = num_gates / num_swap_sim
#     if new_score > current_score:
#         MCTree.nodes[first_selec_node]['sim_score'] = new_score
# =============================================================================

def BackPropagatSimRes(MCTree, start_node, num_swap_sim, num_gates):
    # first way for sim_score
    #h_score = num_gates / num_swap_sim
    #h_score = h_score * 0.7
    # second way for sim_score
    #print(num_gates, num_swap_sim)
    h_score = num_gates * np.power(0.7, num_swap_sim/2)
    
    #print(start_node, num_gates, num_swap_sim, h_score)
    
    # we backpropagation heuristic score
    #MCTree.nodes[start_node]['score'] += h_score
    new_value = MCTree.nodes[start_node]['score'] + h_score
    if new_value > MCTree.nodes[start_node]['global_score']:
        MCTree.nodes[start_node]['global_score'] = new_value
        MCTree.BackPropagationSingleValue(start_node,
                                          value=new_value,
                                          name=None,
                                          mode_BP=['globalscore',[0.7]])
    