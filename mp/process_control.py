#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 14:23:49 2020

@author: xianzhou
"""
import multiprocessing as mp
from circuittransform.method.MCT_Sim_fix_swap import MCTreeSearchSim2
import time

def IniMultiProcessForSim(MCTree,
                          process_num):
    AG = MCTree.AG
    shortest_length_AG = MCTree.shortest_length_AG
    DG = MCTree.DG
    swap_combination = MCTree.swap_combination
    queue_in = mp.Queue()
    queue_out = mp.Queue()
    process_list = []
    sim_con_list = []
    for i in range(process_num):
        flag_con = mp.Value('i', 1)
        #print('start to generate the %dth processor'%i)
        p = mp.Process(target=IniProcessForSim,args=(DG, AG,
                                                     swap_combination,
                                                     shortest_length_AG,
                                                     process_num,
                                                     flag_con,
                                                     queue_in, queue_out,))
        p.start()
        #print('finish to generate the %dth processor'%i)
        process_list.append(p)
        sim_con_list.append(flag_con)
    return process_list, sim_con_list, queue_in, queue_out

def IniProcessForSim(DG, AG, possible_swap_combination,
                     shortest_length_AG,
                     process_num,
                     flag_con,
                     queue_in, queue_out):
    while True:
        try:
            MCT_info = queue_in.get(block=False)
        except:
            pass
        else:
            #t_s = time.time()
            start_nood = MCT_info[0]
            mapping = MCT_info[1]
            selec_num = MCT_info[2]
            root_executed_nodes = MCT_info[3]
            root_executable_nodes = MCT_info[4]
            #print('start sim, node', start_nood)
            #t_s = time.time()
            #print('sim time', selec_num)
            h_score = MCTreeSearchSim2(AG, mapping, DG,
                                       possible_swap_combination,
                                       shortest_length_AG,
                                       root_executable_nodes,
                                       root_executed_nodes,
                                       selec_num)
            #print('finish sim, node', start_nood)
            queue_out.put((start_nood, h_score, selec_num))
            #print('time for sim is ', time.time()-t_s)
            #print('score is', h_score)
        
def StartNewSim(queue_in,
                start_node,
                MCTree,
                selec_num):
    root_executable_nodes = MCTree.nodes[start_node]['executable_vertex']
    root_executed_nodes = MCTree.nodes[start_node]['executed_vertex']
    mapping = MCTree.nodes[start_node]['mapping'].Copy()
    queue_in.put((start_node, mapping, selec_num,
                  root_executed_nodes, root_executable_nodes))
    
# =============================================================================
#     t_s = time.time()
#     h_score = MCTreeSearchSim2(MCTree.AG, mapping, MCTree.DG,
#                                MCTree.swap_combination,
#                                MCTree.shortest_length_AG,
#                                root_executable_nodes.copy(),
#                                root_executed_nodes.copy(),
#                                20)    
#     print('time for sim is ', time.time()-t_s)
#     print('score is', h_score)
#     print('')
# =============================================================================
    
def GetSimRes(MCTree, queue_out,
              finish_all=False,
              sim_count = 1):
    '''obtain sim result and BP the score'''
    #if finish_all == True and sim_count != 0:
    #    print('selec finished but still %d sim remains' %sim_count)
    finish_sim_count = 0
    while sim_count >= 1:
        try:
            start_nood, h_score, selec_num = queue_out.get(block=finish_all)
        except:
            break
        else:
            finish_sim_count += 1
            BackPropagatSimRes(MCTree, h_score, start_nood, selec_num)
            if finish_all == True:
                sim_count -= 1
    return finish_sim_count

def HaltAllSim(p_list, sim_con_list):
    for con in sim_con_list:
        con.value = 0
    for p in p_list:
        #p.join()
        p.terminate()
        
def BackPropagatSimRes(MCTree, h_score, expand_node, selec_num):
    # we backpropagation heuristic score
    g_score_current = MCTree.nodes[expand_node]['global_score']
    #add_visit = sim_swaps / 3 #recommend 3 for fix swaps, 1 for MCTS_sim
    #print('sim node has %d selec' % MCTree.nodes[expand_node]['visited_time_total'])
    #add_visit = max(selec_num-MCTree.nodes[expand_node]['visited_time_total'],
    #                0) * 0.8
    add_visit = 0
    MCTree.nodes[expand_node]['visited_time_total'] += add_visit
    g_score_new = MCTree.nodes[expand_node]['score'] +\
                  h_score*0.7
    current_node = expand_node
    t = 0
    while g_score_new > g_score_current and current_node != MCTree.root_node:
        t += 1
        add_visit = add_visit * 0.2
        '''go to next node'''
        current_node = MCTree.node[current_node]['father_node']
        # update visit count
        MCTree.nodes[current_node]['visited_time_total'] += add_visit
        '''calculate decay rate for propagated score'''          
        g_score_current = MCTree.node[current_node]['global_score']
        g_score_new = g_score_new * 0.7
        g_score_new += MCTree.node[current_node]['score']
        #print('old value', g_score_current)
        #print('new value', g_score_new)
    #print('BP times', t)
    #print('')
# =============================================================================
#     if t == 0:
#         print('new score', g_score_new)
#         print('old scoe', g_score_current)
#         print('old visit', MCTree.nodes[expand_node]['visited_time_total'])
#         print('')
# =============================================================================
    '''renew visited time til root node, this part can be annotated'''
# =============================================================================
#     while current_node != MCTree.root_node:
#         MCTree.node[current_node]['visited_time_total'] += add_visit
#         add_visit = add_visit * 0.7
#         '''go to next node'''
#         current_node = MCTree.node[current_node]['father_node']
# =============================================================================