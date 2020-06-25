# -*- coding: utf-8 -*-
"""
Created on Mon May  4 01:19:27 2020

@author: zxz58
"""
import sys
sys.path.append('..')
from method.remotecnotandwindowbreadth import RemoteCNOTandWindowLookAhead
from circuittransform.operation import GetSubDG
SAHS = RemoteCNOTandWindowLookAhead

def SimSAHS(MCTree, sim_node, sim_cx_num):
    sub_DG = GetSubDG(MCTree.DG,
                      MCTree.nodes[sim_node]['executable_vertex'],
                      MCTree.nodes[sim_node]['executed_vertex'],
                      sim_cx_num)
    res = SAHS(q_phy=None,
               cir_phy=None,
               G=MCTree.AG,
               DG=sub_DG,
               initial_map=MCTree.nodes[sim_node]['mapping'].Copy(),
               shortest_length_Gs=[MCTree.shortest_length_AG, None],
               shortest_path_G=MCTree.shortest_path_AG,
               depth_lookahead=1,
               use_prune=True,
               draw=False,
               DiG=None,
               level_lookahead=None)
    return res[0], len(sub_DG.nodes())