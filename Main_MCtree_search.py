# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 23:21:22 2020

@author: Xiangzhen Zhou

Circuit Transformation via Monte Carlo Tree Search
"""
#############################################################################
'''The following parameters can be adjusted to meet different requirement'''

'''choose a circuits set for benchmark'''
''''small', large'''
benchmark = 'MCTS_small' 

'''repeated times for each circuit'''
repeat_time = 5

'''method for initial mapping'''
'''0: naive; 1: SAHS'''
initial_mapping_control = 1

'''
should we output .qasm files for the transformed circuit?
if True, the output files will be store at ./output_qasm/
'''
output_qasm_file = False

##############################################################################
from inputs.inputgenerator import GenerateArchitectureGraph, CreateDGfromQASMfile
from inputs.shortestpath import ShortestPath
from inputs.map import Map
import networkx as nx
from networkx import DiGraph, Graph
import numpy as np
from qiskit import QuantumCircuit
from qiskit import QuantumRegister
import time
from method.MCtree import MCTree
import matplotlib.pyplot as plt
from mp.process_control_ave_sing_processor import StartNewSim, GetSimRes
from inputs.get_benchmarks import GetBenchmarks, GetMappings

##############################################################################
'''parameter control'''
play_times = [20]
'''mode for search'''
mode = ['fix play', None]
'''mode for simulation'''
mode_sim = [500, 30]

'''initialize parameters'''
# description of architecture graph
num_vertex = 16
# architecture graph generation control
method_AG = ['IBM QX20']


QASM_file = GetBenchmarks(benchmark)
if initial_mapping_control == 1: initial_map_lists = GetMappings(benchmark)


'''generate architecture graph'''
'''q0 - v0, q1 - v1, ...'''
G = GenerateArchitectureGraph(num_vertex, method_AG)
DiG = None
if isinstance(G, DiGraph): #check whether it is a directed graph
    DiG = G
    G = nx.Graph(DiG)
'''calculate shortest path and its length'''
res = ShortestPath(G)
shortest_path_AG = res[1]
shortest_length_AG = (res[0], res[2])
shortest_length_AG = shortest_length_AG[0]

'''only use single swap'''      
possible_swap_combination = []
edges = list(G.edges()).copy()
for current_edge in edges:
    possible_swap_combination.append([current_edge]) 

num_file = 0
compiled_file_name = []

for i in range(len(QASM_file)):
#for i in [0]: #11, 95
    file = QASM_file[i]
    if file[-5:] != '.qasm': file += '.qasm'
    num_file += 1
    '''initial mapping'''
    if initial_mapping_control == 0:
        '''initialize map from logical qubits to physical qubits'''
        '''1-1, 2-2 ...'''
        initial_map_list = list(range(num_vertex))
    if initial_mapping_control == 1:
        initial_map_list = initial_map_lists[i]
    
    print('Number of circuit is', num_file)
    print('Name of circuit is', file)
    for play_time in play_times: #this is selection times for each root node recursion
        mode[1] = play_time
        
        '''MC tree search'''
        '''
        mode: 
            ['fix play', play_num]
            ['entropy play', max_play_num]
            ['fix play with sim', play_num]
        '''
        '''only use single swap'''
                    
        if True:
            res_qasm = CreateDGfromQASMfile(file)
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
        
        for repeat in range(repeat_time):
            print('---------------------------')
            print('The repeated time is ', repeat)
            
            '''initialize map from logical qubits to physical qubits'''
            '''annealing search'''
            initial_map = Map(q_log, G, initial_map_list)
    
            '''MC tree search'''                                                                                      
            tree_size = [1]  
            '''initialize search tree'''
            search_tree = MCTree(G, DG, shortest_length_AG, shortest_path_AG,
                                 possible_swap_combination,
                                 mode_sim)
            s_t1 = time.time()
            search_tree.AddNode(father_node=None, added_SWAP=[], map_new=initial_map)
            search_tree.ExpandNode(0)
            
    
            if mode[0] == 'fix play':
            #    j = -1
                play_num_list = []
                selection_times = mode[1]
                while search_tree.nodes[search_tree.root_node]['num_remain_vertex'] != 0:
            #            j += 1
                    search_depth_c = []
                    '''store pool results for multiprocessing'''
                    sim_count = 0 # current running siumlations
                    sim_count_total = 0# total siumlations
                    selec_times = 0
                    search_tree.simulated_nodes = []
                    #search_tree.CalNodesInDGForSim(search_tree.root_node, mode_sim[1])
                    while selec_times <= selection_times:
                        selec_times += 1
                        '''for each current node, we play rollout for i times'''
                        t_start = time.time()
                        added_depth, expand_node, _, first_selec_node = search_tree.Selection('global_score')
                        search_depth_c.append(added_depth)
                        '''do simulation'''
                        t_start = time.time()
                        if len(search_tree.nodes[expand_node]['son_nodes']) != 0:
                           sim_nodes = [expand_node]
                           for sim_node in sim_nodes:
                                search_tree.CalNodesInDGForSim(sim_node, mode_sim[1])
                                sim_time = mode_sim[0]
                                success = StartNewSim(None,
                                                      sim_node,
                                                      None,
                                                      None,
                                                      search_tree,
                                                      mode_sim[0])
                    '''decision'''
                    search_tree.Decision()
                    #tree_size.append(search_tree.node_size)
                    play_num_list.append(selec_times)
                s_t2 = time.time()
                if len(play_num_list) == 0: play_num_list = [0]
                '''output QASM file'''
                if output_qasm_file == 1:
                    from method.MCTS_draw.MCTS_draw import MCTSToQiskitCir
                    from Qiskitconverter.qasm import CirToQASM
                    from post_processing.qiskit_optimization import PostOptimize
                    cir = MCTSToQiskitCir(search_tree, file, convert_final_map=False)
                    #cir_post = PostOptimize(cir)
                    #print('\nNumber of gates after post optimization is %d' %cir_post.size())
                    #print(cir_post.draw())
                    CirToQASM(cir, file[0:-5]+'_'+str(repeat))
            
            '''display result'''
            total_gate_num = (search_tree.nodes[search_tree.root_node]['num_total_add_gates']*3
                              + cir_log.size())
            print('')
            print('Number of added SWAPs is', search_tree.nodes[search_tree.root_node]['num_total_add_gates'])
            print('Number of added gates is', search_tree.nodes[search_tree.root_node]['num_total_add_gates']*3)
            print('Number of all gates is', total_gate_num)
            print('Total time is', s_t2 - s_t1)
        
    print('===================================')
   