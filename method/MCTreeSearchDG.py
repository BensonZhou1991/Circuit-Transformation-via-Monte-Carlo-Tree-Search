# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 21:31:02 2019

@author: zxz58
"""
import circuittransform as ct
from circuittransform import OperationU, OperationCNOT, OperationSWAP, Map
import time
from circuittransform.method.MCtreeCF import MCTree
from qiskit import QuantumCircuit
from qiskit import QuantumRegister
import matplotlib.pyplot as plt
from circuittransform.map import InitialMapSimulatedAnnealingWeighted
WeightedSimulatedAnnealing = InitialMapSimulatedAnnealingWeighted

num_vertex = 6
display_state = 0 #whether we print the internal states during search process
draw_depth = 0
level_lookahead = [1, 0.8, 0.6, 0.4] #default value [1, 0.8, 0.6, 0.4]

def MCTreeSearchDG(file_name, AG, repeat_time, results, mode, mode_sim,
                     shortest_length_AG, shortest_path_AG,
                     possible_swap_combination,
                     DG=None,
                     initial_map_list=None, draw_circuit=False):
    '''
    mode: 
        ['fix play', play_num]
        ['entropy play', max_play_num]
        ['fix play with sim', play_num]
    '''
    use_SA = False if isinstance(initial_map_list, list) else True
    '''only use single swap'''
    G = AG
    # we comment the following because we decide to input it from outside
# =============================================================================
#     possible_swap_combination = []
#     edges = list(G.edges()).copy()
#     for current_edge in edges:
#         possible_swap_combination.append([current_edge]) 
# =============================================================================
    if DG == None:
        file = file_name
        res_qasm = ct.CreateDGfromQASMfile(file)
        name = file[0:-5]
        
        '''initialize logical quantum circuit'''
        q_log = res_qasm[1][2]
        cir_log = res_qasm[0]
        
        '''initialize physical quantum circuit'''
        q_phy = QuantumRegister(num_vertex, 'v')
        cir_phy = QuantumCircuit(q_phy)
        
        '''generate CNOT operation'''
        total_CNOT = res_qasm[1][3]
     
        '''generate dependency graph'''
        DG = res_qasm[1][0]
    else:
        '''initialize logical quantum circuit'''
        q_log = QuantumRegister(num_vertex, 'q')
        cir_log = QuantumCircuit(q_log)  
        
        '''initialize physical quantum circuit'''
        q_phy = QuantumRegister(num_vertex, 'v')
        cir_phy = QuantumCircuit(q_phy)        

    
#    print('Name of circuit is', file)
    for repeat in range(repeat_time):
        
        '''initialize map from logical qubits to physical qubits'''
        '''annealing search'''
        if use_SA == False:
            t_s = time.time()
            initial_map = Map(q_log, G, initial_map_list)
            t_e = time.time()
        else:
            t_s = time.time()
            map_res = WeightedSimulatedAnnealing(DG,
                                                 AG,
                                                 q_log,
                                                 shortest_length_AG,
                                                 start_map=None,
                                                 convergence=False,
                                                 display=False)
            t_e = time.time()
            initial_map = map_res[0]
            initial_map_list = map_res[1]
        #print('time for initial mapping is', t_e-t_s)

        '''MC tree search'''                                                                                      
        tree_size = [1]
        search_depth_ave = []
        search_depth_max = []
        search_depth_min = []
        names_father = ['executable_vertex', 'executed_vertex']
        names_son = ('added_SWAP', 'visited_time', 'visited_time_total',
                     'score', 'global_score')
        s_t1 = time.time()
        '''initialize search tree'''
        search_tree = MCTree(G, DG, shortest_length_AG, shortest_path_AG,
                             possible_swap_combination,
                             level_lookahead,
                             mode_sim)
        search_tree.AddNode(father_node=None, added_SWAP=[], map_new=initial_map)
        search_tree.ExpandNode(0)
        
    # =============================================================================
    #         for i in range(200):
    #             '''for each current node, we play rollout for i times'''
    #             search_depth.append(search_tree.Selection('global_score'))
    #         search_tree.PrintSonNodesArgs(0, names_son)
    #         search_tree.Decision('global_score')
    #         tree_size.append(search_tree.node_size)
    # =============================================================================

        if mode[0] == 'fix play':
        #    j = -1
            #results[name]['ave play num'].append(mode[1])
            while search_tree.node[search_tree.root_node]['num_remain_vertex'] != 0:
        #            j += 1
                search_depth_c = []
                for i in range(mode[1]):
                    '''for each current node, we play rollout for i times'''
                    added_pepth = search_tree.Selection('global_score')
                    search_depth_c.append(added_pepth)
        #            if j == 10: search_tree.PrintSonNodesArgs(search_tree.root_node, names_son)
                if display_state == 1:
                    search_tree.PrintSonNodesArgs(search_tree.root_node, names_son)
                search_tree.Decision('global_score')
                search_depth_ave.append(sum(search_depth_c)/len(search_depth_c))
                search_depth_max.append(max(search_depth_c))
                search_depth_min.append(min(search_depth_c))
                #tree_size.append(search_tree.node_size)
            s_t2 = time.time()

        if mode[0] == 'fix play with sim':
        #    j = -1
            results[name]['ave play num'].append(mode[1])
            while search_tree.node[search_tree.root_node]['num_remain_vertex'] != 0:
        #            j += 1
                search_depth_c = []
                for i in range(mode[1]):
                    '''for each current node, we play rollout for i times'''
                    search_depth_c.append(search_tree.SelectionWithSimulation('global_score'))
        #            if j == 10: search_tree.PrintSonNodesArgs(search_tree.root_node, names_son)
                    
                search_tree.Decision('global_score')
                search_depth_ave.append(sum(search_depth_c)/len(search_depth_c))
                search_depth_max.append(max(search_depth_c))
                search_depth_min.append(min(search_depth_c))
                #tree_size.append(search_tree.node_size)
            s_t2 = time.time()
        
        if mode[0] == 'entropy play':
            '''use entropy'''
            play_time = []
            while search_tree.node[search_tree.root_node]['num_remain_vertex'] != 0:
    #            j += 1
                search_depth_c = []
                k = search_tree.CalUnConfidence(search_tree.root_node, 'global_score')
                #print('%.2f'%k)
                k = int(mode[1] * k)
                k = max(k, 10)
                play_time.append(k)
                for i in range(k):
                    '''for each current node, we play rollout for i times'''
                    search_depth_c.append(search_tree.Selection('global_score'))
    #            if j == 10: search_tree.PrintSonNodesArgs(search_tree.root_node, names_son)
                search_tree.Decision('global_score')
                search_depth_ave.append(sum(search_depth_c)/len(search_depth_c))
                search_depth_max.append(max(search_depth_c))
                search_depth_min.append(min(search_depth_c))
                #tree_size.append(search_tree.node_size)
            s_t2 = time.time()
            ave_play_time = sum(play_time) / (len(play_time)+0.00001)
            results[name]['ave play num'].append(ave_play_time)
            #print('Average play time is', ave_play_time)
        
        '''display result'''
        total_gate_num = (search_tree.node[search_tree.root_node]['num_total_add_gates']*3
                          + cir_log.size())
        add_gate_num = search_tree.node[search_tree.root_node]['num_total_add_gates']*3
        final_map_list = search_tree.node[search_tree.root_node]['mapping']
        #print('Number of total gate number is', total_gate_num)
        return add_gate_num, initial_map_list, final_map_list.MapToList()
# =============================================================================
#         print('Number of added SWAPs is', search_tree.node[search_tree.root_node]['num_total_add_gates'])
#         print('Number of added gates is', search_tree.node[search_tree.root_node]['num_total_add_gates']*3)
#         print('Number of total gate number is', total_gate_num)
#         print('First finished added gates is', search_tree.first_finish_add_swaps*3)
#         print('First finished time is', search_tree.first_finish_time - s_t1)
#         print('Total time is', s_t2 - s_t1)
#         print('Size of tree is', search_tree.node_count)
# =============================================================================
