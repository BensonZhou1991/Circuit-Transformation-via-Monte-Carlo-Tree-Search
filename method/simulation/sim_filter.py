# -*- coding: utf-8 -*-
"""
Created on Mon May  4 02:41:42 2020

@author: zxz58
"""
import networkx as nx
from circuittransform.Li.QTokyo_f import qct

def SimFilter(MCTree, sim_node):
    sim_nodes = MCTree.nodes[sim_node]['sim_nodes']
    num_cx = len(sim_nodes)
    C = []
    EG = nx.edges(MCTree.AG)
    for node in sim_nodes:
        C.append(MCTree.DG.nodes[node]['operation'].involve_qubits_list)
    initial_map=MCTree.nodes[sim_node]['mapping']
    #print(C)
    #print(initial_map.MapToList())
    tau = initial_map.MapToListReverse()
    #print(tau)
    #print()
    C_out, cost_time = qct(C, MCTree.AG, EG, tau, '12')
    num_swap = (len(C_out) - num_cx) / 3
    return num_swap, num_cx