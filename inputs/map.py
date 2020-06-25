# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 22:22:20 2019

@author: zxz58

This module is for the map from logical qubits to physical qubits
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

class Map(object):
    
    instances_count = 0
    
    def __init__(self, q_reg, G, initial_map = None):
        '''
        if initial_map is None, q[0] -> v[0], q[1] -> v[1]...
        else, initial_map is a list q[0] -> initial_map[0], q[1] -> initial_map[1]...
        
        Input:
            q_reg: logical qubits represented by QuantumRegister in Qiskit
            initial_map: a dic{qi: vj...} reprenenting the initial map
                        or a list [v0, v1...] whose entry reprenents the physical qubit
        
        _DomToCod: dic from domain of definition to codomain
        _CodToDom: dic from codomain to domain of definition
        '''
        vertex = list(G.nodes())
        num_q = len(q_reg)
        self.num_v = len(list(G))
        self.logical_quantum_register = q_reg
        self.architecture_graph = G
        self.__DomToCod = {}
        self.__CodToDom = {}
        if initial_map == None:
            for q in q_reg:
                self.__DomToCod[q] = vertex[q[1]]
                self.__CodToDom[vertex[q[1]]] = q
            # define empty vertexes in architecture graph
            if self.num_v > num_q:
                for i in range(num_q, self.num_v):
                    self.__CodToDom[vertex[i]] = []
        else:
            if isinstance(initial_map, dict):
                # initialize __CodToDom
                for v in vertex:
                    self.__CodToDom[v] = []
                keys = list(initial_map.keys())
                for q in keys:
                    self.__DomToCod[q] = initial_map[q]
                    self.__CodToDom[initial_map[q]] = q
            if isinstance(initial_map, list):
                # initialize __CodToDom
                for v in vertex:
                    self.__CodToDom[v] = []
                for i in range(len(q_reg)):
                    q = q_reg[i]
                    self.__DomToCod[q] = initial_map[i]
                    self.__CodToDom[initial_map[i]] = q

        Map.instances_count = Map.instances_count + 1 
    
    def Copy(self):
        '''return a copy of a Map instance'''
        return Map(self.logical_quantum_register, self.architecture_graph, self.__DomToCod)
        
    def DomToCod(self, qubit):
        if isinstance(qubit, int):
            qubit = self.logical_quantum_register[qubit]
        else:
            qubit = self.logical_quantum_register[qubit.index]
        return self.__DomToCod[qubit]
    
    def CodToDom(self, vertex):
        return self.__CodToDom[vertex]

    def LogToPhy(self, qubit):
        if isinstance(qubit, int):
            qubit = self.logical_quantum_register[qubit]
        else:
            qubit = self.logical_quantum_register[qubit.index]
        return self.__DomToCod[qubit]
    
    def PhyToLog(self, vertex):
        return self.__CodToDom[vertex]
    
    def RenewMapViaExchangeDom(self, q0, q1):
        '''
        exchange mapping of qubits q0 and q1
        '''
        temp_q0 = self.__DomToCod[q0]
        temp_q1 = self.__DomToCod[q1]
        self.__DomToCod[q0] = temp_q1
        self.__DomToCod[q1] = temp_q0
        self.__CodToDom[temp_q0] = q1
        self.__CodToDom[temp_q1] = q0

    def RenewMapViaExchangeCod(self, v0, v1):
        '''
        exchange mapping of vertexes (physical qubits) v0 and v1
        '''
        temp_v0 = self.__CodToDom[v0]
        temp_v1 = self.__CodToDom[v1]
        self.__CodToDom[v0] = temp_v1
        self.__CodToDom[v1] = temp_v0
        if temp_v0 != []:
            self.__DomToCod[temp_v0] = v1
        if temp_v1 != []:
            self.__DomToCod[temp_v1] = v0
    
    def RenewMapViaSWAP(self, operation_swap):
        '''
        update map after swap operation
        '''
        q0 = operation_swap.involve_qubits[0]
        q1 = operation_swap.involve_qubits[1]
        self.RenewMapViaExchangeDom(q0, q1)
        
    def MapToTuple(self):
        '''
        calculate corresponding tuple for the map
        it can be used to indetify the identical mapping
        '''
        out = []
        for q in self.logical_quantum_register:
            out.append(self.DomToCod(q))
        
        return out
    
    def MapToList(self):
        '''
        return the list (index is log) of the map
        '''
        map_list = []
        q_log = self.logical_quantum_register
        for q in q_log:
            map_list.append(self.LogToPhy(q))
        
        return map_list
    
    def MapToListReverse(self):
        '''
        return the list (index is phy) of the map
        '''
        map_list = []
        for q_phy in range(self.num_v):
            q_log = self.PhyToLog(q_phy)
            if q_log == []: q_log = 20
            if not isinstance(q_log, int):
                q_log = q_log.index
            map_list.append(q_log)
        
        return map_list        
    
def SortKey(qubit_weight):
    return qubit_weight[1]

def FindInitialMapping(DG, q_log, G, shortest_length_G):
    '''
    find a initial mapping
    step1: count the used times for each logical qubit and logical qubits pairs
    step2: choose the logical qubit mostly used and allocate it to physical qubit with the best connectivity
    step3: choose the logical qubits who will be used mostly to compose to CNOT with qubit in step2 and allocate
    them to to physical qubits conncted to the one in step2
    input:
        G must be an undirected graph
    '''
    len_q = len(q_log)
    len_v = len(G.nodes())
    qubit_weight = []
    vertex_weight = []
    qubit_qubit_weight = []
    initial_map = [-1]*len_q
    for q in range(len_q):
        qubit_weight.append([q, 0])
        attend_list = []
        for q2 in range(len_q):
            attend_list.append([q2, 0])
        qubit_qubit_weight.append(attend_list)
    for v in range(len_v):
        vertex_weight.append([v, 0])
    '''calculate the used times for each logical qubit'''
    for node in DG.nodes():
        op = DG.nodes[node]['operation']
        involve_qubits = op.involve_qubits
        qubit_qubit_weight[involve_qubits[0][1]][involve_qubits[1][1]][1] += 1
        qubit_qubit_weight[involve_qubits[1][1]][involve_qubits[0][1]][1] += 1
        for q_used in involve_qubits:
            qubit_weight[q_used[1]][1] += 1
    '''calculate the performance for each physical qubit'''
    for node in G.nodes():
        vertex_weight[node][1] += sum(shortest_length_G[node].values())
    '''sort the logical qubits'''
    #print(qubit_weight)
    qubit_weight.sort(key=SortKey, reverse=True)
    vertex_weight.sort(key=SortKey, reverse=False)
    
    #print(qubit_weight)
    #print(qubit_qubit_weight)
    #print(vertex_weight)
    
    for [current_qubit, current_qubit_weight] in qubit_weight:
        if initial_map[current_qubit] != -1: continue
        '''search best node in AG'''
        for [best_vertex, best_vertex_weight] in vertex_weight:
            if not best_vertex in initial_map: break
        '''allocate this node to qubit'''
        initial_map[current_qubit] = best_vertex
        '''allocate nodes in AG to qubits adjacent to current one'''
        next_qubits_list = qubit_qubit_weight[current_qubit]
        next_qubits_list.sort(key=SortKey, reverse=True)
        for [next_qubit, next_weight] in next_qubits_list:
            if next_weight == 0: break # if this qubit has no connection to current one
            if initial_map[next_qubit] != -1: continue
            next_vertexes = list(G[best_vertex])
            for [best_next_vertex, best_next_vertex_weight] in vertex_weight:
                if not best_next_vertex in next_vertexes: continue
                if not best_next_vertex in initial_map:
                    initial_map[next_qubit]
    return Map(q_log, G, initial_map), initial_map

def initpara():
    '''only for simulated annealing'''
    alpha = 0.98
    t = (1,100)#(1,100)
    markovlen = 100
    return alpha,t,markovlen

def CalCost(map_list, DG, counted_CNOT_nodes, shortest_length_G, shortest_path_G, q_log, G, DiG):
    '''only for simulated annealing'''
    print('Deprecated function!')
    if isinstance(map_list, list):
        cost = ct.HeuristicCostZulehner(Map(q_log, G, map_list), DG, counted_CNOT_nodes, shortest_length_G, shortest_path_G, DiG)
    else:
        cost = ct.HeuristicCostZulehner(map_list, DG, counted_CNOT_nodes, shortest_length_G, shortest_path_G, DiG)
    return cost[1]

def MapListReverse(map_list, num_v):
    v_list = [-1] * num_v
    for i in range(len(map_list)):
        v_list[map_list[i]] = i
    return v_list

def VListReverse(v_list, num_q):
    map_list = [-1] * num_q
    for i in range(len(v_list)):
        if v_list[i] == -1: continue
        map_list[v_list[i]] = i
    return map_list

def InitialMapSimulatedAnnealing(start_map, DG, G, DiG, q_log, shortest_length_G, shortest_path_G, num_consider_gates=0, convergence=False):
    '''
    this function is modified from "https://blog.csdn.net/qq_34798326/article/details/79013338"
    '''
    AG = G
    if convergence == True:
        temp = []
        solution = []
        solution_best = []
    if len(start_map) != len(G.nodes()):
        for v in G.nodes():
            if not v in start_map: start_map.append(v)
    if num_consider_gates <= 1:
        num_consider_gates = len(DG.nodes()) * num_consider_gates
        if num_consider_gates < 50: num_consider_gates = 50
    if num_consider_gates > len(DG.nodes()): num_consider_gates = len(DG.nodes())
    DG_copy = DG.copy()
    '''generate counted CNOTs'''
    counted_CNOT_nodes = []
    removed_nodes = list(range(len(DG)))
    executable_vertex = ct.FindExecutableNode(DG)
    executed_vertex = []
    while len(counted_CNOT_nodes) < num_consider_gates:        
        counted_CNOT_nodes.extend(executable_vertex)
        for r_node in executable_vertex: removed_nodes.remove(r_node)
        # go to next level
        executable_vertex = ct.FindExecutableNode(DG,
                                                  executed_vertex,
                                                  executable_vertex,
                                                  executable_vertex.copy())   
    DG_copy.remove_nodes_from(removed_nodes)
    print('number of counted gates is', len(counted_CNOT_nodes))
    '''gen cost matrix'''
    num_q = len(AG.nodes())
    cost_m = InitialCostMatrixWeighted(DG_copy, num_q, method=4)
    '''Simulated Annealing'''
    solutionnew = start_map
    num = len(start_map)
    #valuenew = np.max(num)
    solutioncurrent = solutionnew.copy()
    valuecurrent = 99000  #np.max这样的源代码可能同样是因为版本问题被当做函数不能正确使用，应取一个较大值作为初始值
    
    #print(valuecurrent)
    
    solutionbest = solutionnew.copy()
    valuebest = 99000 #np.max
    alpha,t2,markovlen = initpara()
    t = t2[1]#temperature
    result = [] #记录迭代过程中的最优解
    
    while t > t2[0]:
        for i in np.arange(markovlen):
            #下面的两交换和三角换是两种扰动方式，用于产生新解
            if np.random.rand() > 0.5:# 交换路径中的这2个节点的顺序
                # np.random.rand()产生[0, 1)区间的均匀随机数
                while True:#产生两个不同的随机数
                    loc1 = np.int(np.around(np.random.rand()*(num-1)))
                    loc2 = np.int(np.around(np.random.rand()*(num-1)))
                    ## print(loc1,loc2)
                    if loc1 != loc2:
                        break
                solutionnew[loc1],solutionnew[loc2] = solutionnew[loc2],solutionnew[loc1]
            else: #三交换
                while True:
                    loc1 = np.int(np.around(np.random.rand()*(num-1)))
                    loc2 = np.int(np.around(np.random.rand()*(num-1))) 
                    loc3 = np.int(np.around(np.random.rand()*(num-1)))
                    if((loc1 != loc2)&(loc2 != loc3)&(loc1 != loc3)):
                        break
                # 下面的三个判断语句使得loc1<loc2<loc3
                if loc1 > loc2:
                    loc1,loc2 = loc2,loc1
                if loc2 > loc3:
                    loc2,loc3 = loc3,loc2
                if loc1 > loc2:
                    loc1,loc2 = loc2,loc1
                #下面的三行代码将[loc1,loc2)区间的数据插入到loc3之后
                tmplist = solutionnew[loc1:loc2].copy()
                solutionnew[loc1:loc3-loc2+1+loc1] = solutionnew[loc2:loc3+1].copy()
                solutionnew[loc3-loc2+1+loc1:loc3+1] = tmplist.copy()  
            valuenew = CalCostMatrixWeighted(cost_m, solutionnew,
                                             shortest_length_G, num_q)
           # print (valuenew)
            if valuenew < valuecurrent: #接受该解
                #更新solutioncurrent 和solutionbest
                valuecurrent = valuenew
                solutioncurrent = solutionnew.copy()
                #renew best solution
                if valuenew < valuebest:
                    valuebest = valuenew
                    solutionbest = solutionnew.copy()
            else:#按一定的概率接受该解
                if np.random.rand() < np.exp(-(valuenew-valuecurrent)/t):
                    valuecurrent = valuenew
                    solutioncurrent = solutionnew.copy()
                else:
                    solutionnew = solutioncurrent.copy()

            if convergence == True:
                temp.append(t)
                solution.append(valuecurrent)
                solution_best.append(valuebest)
        
        t = alpha*t
        #print(valuebest)
        result.append(valuebest)
    if display == True:
        print('initial_map is', solutionbest)
        print('best value for SA is', valuebest)
    '''draw convergence graph'''
    if convergence == True:
        figure_fig = plt.figure()
        plt.grid()
        plt.xlabel('Times of Iteration')
        plt.ylabel('Cost of States')
        plt.plot(solution)
        plt.plot(solution_best)
        figure_fig.savefig('simulated annealing convergence.eps', format='eps', dpi=1000)
        
    return Map(q_log, G, solutionbest), solutionbest

def InitialCostMatrixWeighted(DG, num_q, method=2, add_weight=False):
    cost_m = np.zeros((num_q, num_q))
    if method == 1:
        '''method 1, exponential decay'''
        weight = 1
        decay = 0.99 # default 0.99
        while len(DG.nodes()) != 0:
            current_nodes = ct.FindExecutableNode(DG)
            for node in current_nodes:
                op = DG.nodes[node]['operation']
                qubits = op.involve_qubits
                cost_m[qubits[0][1]][qubits[1][1]] += weight
                #cost_m[qubits[1][1]][qubits[0][1]] += weight
            DG.remove_nodes_from(current_nodes)
            weight = weight * decay
    if method == 2:
        '''method 2, linear decay'''
        num_CX = len(DG.nodes())
        num_CX_current = num_CX
        weight = 1
        while len(DG.nodes()) != 0:
            weight = num_CX_current / num_CX
            current_nodes = FindExecutableNode(DG)
            num_CX_current -= len(current_nodes)
            for node in current_nodes:
                op = DG.nodes[node]['operation']
                if add_weight == True:
                    flag = 1
                    '''if comment the following, we ignore the successive CX'''
                    if DG.out_degree(node) == 1:
                        qubits = op.ToTuple()
                        op_next = DG.nodes[list(DG.successors(node))[0]]['operation']
                        qubits_next = op_next.ToTuple()
                        if qubits[0] == qubits_next[0] and qubits[1] == qubits_next[1]:
                            flag = 0
                        if qubits[0] == qubits_next[1] and qubits[1] == qubits_next[0]:
                            flag = 0 
                        
                    op.weight = weight * flag
                qubits = op.involve_qubits
                if isinstance(qubits[0], int):
                    qubits[0] = qubits[0], qubits[0]
                    qubits[1] = qubits[1], qubits[1]
                if add_weight == True:
                    cost_m[qubits[0].index][qubits[1].index] += op.weight
                else:
                    cost_m[qubits[0].index][qubits[1].index] += weight
                #cost_m[qubits[1][1]][qubits[0][1]] += weight
            DG.remove_nodes_from(current_nodes)
    if method == 3:
        '''method 3, weighted cos decay'''
        num_CX = len(DG.nodes())
        num_CX_current = 0
        weight = 1
        while len(DG.nodes()) != 0:
            weightt = num_CX_current / num_CX
            weightt = np.power(weightt, 1)
            weight = (np.cos(np.pi * weightt)+1) / 2
            print(weight)
            current_nodes = ct.FindExecutableNode(DG)
            num_CX_current += len(current_nodes)
            for node in current_nodes:
                op = DG.nodes[node]['operation']
                qubits = op.involve_qubits
                cost_m[qubits[0][1]][qubits[1][1]] += weight
                #cost_m[qubits[1][1]][qubits[0][1]] += weight
            DG.remove_nodes_from(current_nodes)
    if method == 4:
        '''method 4, no decay'''
        num_CX = len(DG.nodes())
        num_CX_current = num_CX
        weight = 1
        while len(DG.nodes()) != 0:
            current_nodes = ct.FindExecutableNode(DG)
            num_CX_current -= len(current_nodes)
            for node in current_nodes:
                op = DG.nodes[node]['operation']
                if add_weight == True:
                    #print(DG.out_degree(node))
                    flag = 1
                    '''if comment the following, we ignore the successive CX'''
                    if DG.out_degree(node) == 1:
                        qubits = op.ToTuple()
                        op_next = DG.nodes[list(DG.successors(node))[0]]['operation']
                        qubits_next = op_next.ToTuple()
                        if qubits[0] == qubits_next[0] and qubits[1] == qubits_next[1]:
                            flag = 0
                        if qubits[0] == qubits_next[1] and qubits[1] == qubits_next[0]:
                            flag = 0                        
                    
                    op.weight = weight * flag
                qubits = op.involve_qubits
                if isinstance(qubits[0], int):
                    qubits[0] = qubits[0], qubits[0]
                    qubits[1] = qubits[1], qubits[1]
                if add_weight == True:
                    cost_m[qubits[0][1]][qubits[1][1]] += op.weight
                else:
                    cost_m[qubits[0][1]][qubits[1][1]] += weight
                #cost_m[qubits[1][1]][qubits[0][1]] += weight
            DG.remove_nodes_from(current_nodes)
    return cost_m

def CalCostMatrixWeighted(cost_m, current_sol, shortest_length_G, num_q):
    cost_total = 0
    for q1_log in range(num_q):
        for q2_log in range(num_q):
            q1_phy, q2_phy = current_sol[q1_log], current_sol[q2_log]
            num_swap = shortest_length_G[q1_phy][q2_phy] - 1
            cost_total += num_swap * cost_m[q1_log][q2_log]
    return cost_total          

def InitialMapSimulatedAnnealingWeighted(DG,
                                         AG,
                                         q_log,
                                         shortest_length_G,
                                         start_map=None,
                                         convergence=False,
                                         display=True):
    '''
    this function is modified from "https://blog.csdn.net/qq_34798326/article/details/79013338"
    '''
    num_q = len(AG.nodes()) # num of physical qubits
    G = AG
    if convergence == True:
        temp = []
        solution = []
        solution_best = []
    if start_map == None: start_map = list(range(num_q))
    if len(start_map) != len(G.nodes()):
        '''
        if logical qubits is less than physical, we extend logical qubit to
        ensure the completeness and delete added qubits at the end of the
        algorithm
        '''
        for v in G.nodes():
            if not v in start_map: start_map.append(v)
    DG_copy = DG.copy()
    '''gen cost matrix'''
    cost_m = InitialCostMatrixWeighted(DG_copy, num_q)
    '''Simulated Annealing'''
    solutionnew = start_map
    num = len(start_map)
    #valuenew = np.max(num)
    solutioncurrent = solutionnew.copy()
    valuecurrent = 99000  #np.max这样的源代码可能同样是因为版本问题被当做函数不能正确使用，应取一个较大值作为初始值
    
    #print(valuecurrent)
    
    solutionbest = solutionnew.copy()
    valuebest = 99000 #np.max
    alpha,t2,markovlen = initpara()
    t = t2[1]#temperature
    result = [] #记录迭代过程中的最优解
    
    while t > t2[0]:
        for i in np.arange(markovlen):
            #下面的两交换和三角换是两种扰动方式，用于产生新解
            if np.random.rand() > 0.5:# 交换路径中的这2个节点的顺序
                # np.random.rand()产生[0, 1)区间的均匀随机数
                while True:#产生两个不同的随机数
                    loc1 = np.int(np.around(np.random.rand()*(num-1)))
                    loc2 = np.int(np.around(np.random.rand()*(num-1)))
                    ## print(loc1,loc2)
                    if loc1 != loc2:
                        break
                solutionnew[loc1],solutionnew[loc2] = solutionnew[loc2],solutionnew[loc1]
            else: #三交换
                while True:
                    loc1 = np.int(np.around(np.random.rand()*(num-1)))
                    loc2 = np.int(np.around(np.random.rand()*(num-1))) 
                    loc3 = np.int(np.around(np.random.rand()*(num-1)))
                    if((loc1 != loc2)&(loc2 != loc3)&(loc1 != loc3)):
                        break
                # 下面的三个判断语句使得loc1<loc2<loc3
                if loc1 > loc2:
                    loc1,loc2 = loc2,loc1
                if loc2 > loc3:
                    loc2,loc3 = loc3,loc2
                if loc1 > loc2:
                    loc1,loc2 = loc2,loc1
                #下面的三行代码将[loc1,loc2)区间的数据插入到loc3之后
                tmplist = solutionnew[loc1:loc2].copy()
                solutionnew[loc1:loc3-loc2+1+loc1] = solutionnew[loc2:loc3+1].copy()
                solutionnew[loc3-loc2+1+loc1:loc3+1] = tmplist.copy()  
            valuenew = CalCostMatrixWeighted(cost_m, solutionnew,
                                             shortest_length_G, num_q)
           # print (valuenew)
            if valuenew < valuecurrent: #接受该解
                #更新solutioncurrent 和solutionbest
                valuecurrent = valuenew
                solutioncurrent = solutionnew.copy()
                #renew best solution
                if valuenew < valuebest:
                    valuebest = valuenew
                    solutionbest = solutionnew.copy()
            else:#按一定的概率接受该解
                if np.random.rand() < np.exp(-(valuenew-valuecurrent)/t):
                    valuecurrent = valuenew
                    solutioncurrent = solutionnew.copy()
                else:
                    solutionnew = solutioncurrent.copy()

            if convergence == True:
                temp.append(t)
                solution.append(valuecurrent)
                solution_best.append(valuebest)
        
        t = alpha*t
        #print(valuebest)
        result.append(valuebest)
    if display == True:
        print('initial_map is', solutionbest)
        print('best value for SA is', valuebest)
    '''draw convergence graph'''
    if convergence == True:
        figure_fig = plt.figure()
        plt.grid()
        plt.xlabel('Times of Iteration')
        plt.ylabel('Cost of States')
        plt.plot(solution)
        plt.plot(solution_best)
        figure_fig.savefig('simulated annealing convergence.eps', format='eps', dpi=1000)
        
    return Map(q_log, G, solutionbest), solutionbest

def MapListPhyToLog(map_list, q_phy_list):
    q_log_list = [-1] * len(q_phy_list)
    if len(q_phy_list) == 2:
        for q_log in range(len(map_list)):
            if map_list[q_log] == q_phy_list[0]:
                q_log_list[0] = q_log
            if map_list[q_log] == q_phy_list[1]:
                q_log_list[1] = q_log
    if len(q_phy_list) == 1:
        for q_log in range(len(map_list)):
            if map_list[q_log] == q_phy_list[0]:
                q_log_list[0] = q_log  
    return q_log_list



def InitialMapAStar(DG, AG, num_q_log, num_q_phy, shortest_length_AG):
    '''only for undirected AG like Q20'''
    ini_map_list = [-1] * num_q_phy
    num_mapped_q = 0
    edges_AG = AG.edges()
    executable_vertex = ct.FindExecutableNode(DG)
    executed_vertex = []
    while num_mapped_q < num_q_log and len(executable_vertex) > 0:
        #print(executable_vertex)
        for node in executable_vertex:
            flag_find_mapping = False
            op = DG.nodes[node]['operation']
            qubits_in_log = op.input_qubits #input logical qubits
            qubits_in_log = qubits_in_log[0][1], qubits_in_log[1][1]
            qubits_in_phy = (ini_map_list[qubits_in_log[0]],
                             ini_map_list[qubits_in_log[1]])
            #print(qubits_in_log)
            #print(qubits_in_phy)
            if qubits_in_phy[0] == -1 and qubits_in_phy[1] == -1:
                for edge in edges_AG:
                    q_log_list = MapListPhyToLog(ini_map_list, edge)
                    if q_log_list == [-1, -1]:
                        ini_map_list[qubits_in_log[0]] = edge[0]
                        ini_map_list[qubits_in_log[1]] = edge[1]
                        flag_find_mapping = True
                        continue
                if flag_find_mapping == False:
                    raise(Exception('Cant find any vacant edge'))
                num_mapped_q += 2
            if qubits_in_phy[0] != -1 and qubits_in_phy[1] == -1:
                #print(qubits_in_phy)
                min_dis = 1000
                chosen_q_phy = None
                for q_phy in range(num_q_phy):
                    q_log = MapListPhyToLog(ini_map_list, [q_phy])[0]
                    if q_log == -1:
                        currant_dis = shortest_length_AG[qubits_in_phy[0]][q_phy]
                        if currant_dis < min_dis:
                            chosen_q_phy = q_phy
                            min_dis = currant_dis
                if chosen_q_phy == None:
                    raise(Exception('Cant find any vacant edge'))
                ini_map_list[qubits_in_log[1]] = chosen_q_phy
                num_mapped_q += 1
                            
            if qubits_in_phy[0] == -1 and qubits_in_phy[1] != -1:
                #print(qubits_in_phy)
                min_dis = 1000
                chosen_q_phy = None
                for q_phy in range(num_q_phy):
                    q_log = MapListPhyToLog(ini_map_list, [q_phy])[0]
                    if q_log == -1:
                        currant_dis = shortest_length_AG[qubits_in_phy[1]][q_phy]
                        if currant_dis < min_dis:
                            chosen_q_phy = q_phy
                            min_dis = currant_dis
                if chosen_q_phy == None:
                    raise(Exception('Cant find any vacant edge'))
                ini_map_list[qubits_in_log[0]] = chosen_q_phy
                num_mapped_q += 1
        ct.FindExecutableNode(DG,
                              executed_vertex,
                              executable_vertex,
                              executable_vertex)
    return ini_map_list

def FindConnectedUnoccupiedNode(AG, map_list, num_q):
    q_candidates = []
    for edge in AG.edges():
        if (not edge[0] in map_list) and (not edge[1] in map_list):
            if not edge[0] in q_candidates:
                q_candidates.append(edge[0])
            if not edge[1] in q_candidates:
                q_candidates.append(edge[1])
    for q_phy in q_candidates:
        for q_phy2 in map_list:
            if (q_phy2, q_phy) in AG.edges:
                return q_phy
    for q_phy in range(num_q):
        if not q_phy in map_list:
            for q_phy2 in map_list:
                if (q_phy2, q_phy) in AG.edges:
                    return q_phy

def MapConvert2(tau, file_name, num_q, shortest_length_AG):
    '''Convert tau to a list where index represents log and value phy'''
    map_list = MapConvert(tau)
    res_qasm = ct.CreateDGfromQASMfile(file_name)
    '''gen max length'''
    shortest_length = shortest_length_AG
    max_length = 0
    for node1 in range(num_q):
        for node2 in range(num_q):
            if shortest_length[node1][node2] > max_length:
                max_length = shortest_length[node1][node2]
    '''generate dependency graph'''
    DG = res_qasm[1][0]
    executable_nodes = ct.FindExecutableNode(dependency_graph=DG,
                                             executed_vertex=None,
                                             executable_vertex=None,
                                             removed_vertexes=None)
    
    executed_nodes = []
    removed_nodes = []
    num_q_log = DG.num_q_log
    '''get allocated logical qubit number'''
    num_q_log_allocated = 0
    for q in range(num_q):
        if tau[q] < num_q: num_q_log_allocated += 1
    while num_q_log_allocated < num_q_log:
        #print(num_q_log_allocated, num_q_log)
        '''execute nodes'''       
        flag_finish = False
        while flag_finish == False:
            flag_finish = True
            #print('executable', executable_nodes)
            #print('removed', removed_nodes)
            executable_nodes = ct.FindExecutableNode(DG,
                                                     executed_nodes,
                                                     executable_nodes,
                                                     removed_nodes)
            removed_nodes = []
            '''find node whose distance is 1'''
            for node in executable_nodes:
                qubits = DG.nodes[node]['operation'].InvolveQubitsList()
                q_phy1, q_phy2 = map_list[qubits[0]], map_list[qubits[1]]
                if q_phy1 != -1 and q_phy2 != -1 and not node in removed_nodes:
                    removed_nodes.append(node)
                    flag_finish = False
            #print('executable', executable_nodes)
        '''allocate one qubit'''
        dis_deduct = -1
        chosen_log_q = -1
        chosen_phy_q = -1
        chosen_node = -1
        
        for q_log in range(num_q_log):
            if map_list[q_log] == -1:
                # this logical qubit is unallocated
                for node in executable_nodes:
                    # find node involves this qubit
                    qubits = DG.nodes[node]['operation'].InvolveQubitsList()
                    if (q_log == qubits[0] and map_list[qubits[1]] != -1) or \
                       (q_log == qubits[1] and map_list[qubits[0]] != -1):
                           # successfully find an involved CNOT
                           # one of input qubits is allocated and the other not
                           if q_log == qubits[0]:
                                q_phy1 = map_list[qubits[1]]
                           else:
                                q_phy1 = map_list[qubits[0]]
                           for q_phy2 in range(num_q):
                               # evaluate different physical qubits
                               if not q_phy2 in map_list:
                                   #this phy q is unallocated
                                   dis_deduct_c = max_length - shortest_length[q_phy2][q_phy1]
                                   if dis_deduct_c > dis_deduct:
                                       dis_deduct = dis_deduct_c
                                       chosen_log_q = q_log
                                       chosen_phy_q = q_phy2
                                       chosen_node = node
                    if (q_log == qubits[0] and map_list[qubits[1]] == -1) or \
                       (q_log == qubits[1] and map_list[qubits[0]] == -1):
                           # both qubits is unallocated
                           chosen_phy_q = FindConnectedUnoccupiedNode(AG, map_list, num_q)
                           dis_deduct = max_length
                           chosen_log_q = q_log
                           chosen_node = None
        #print('chosen p_log', chosen_log_q)
        #print('chosen p_phy', chosen_phy_q)
        map_list[chosen_log_q] = chosen_phy_q
        num_q_log_allocated += 1
        removed_nodes = [chosen_node]
        
    return map_list                            

def MapConvert(tau):
    '''Convert tau to a list where index represents log and value phy'''
    num_q = len(tau)
    map_list = [-1] * num_q
    for q_phy in range(num_q):
        q_log = tau[q_phy]
        #map_list[q_log] = q_phy
        if q_log < num_q: map_list[q_log] = q_phy
    return map_list

def FinishMap(map_list):
    num_q = len(map_list)
    for q_log in range(num_q):
        q_phy = map_list[q_log]
        if q_phy == -1:
            for q_phy in range(num_q):
                if not q_phy in map_list:
                    map_list[q_log] = q_phy
    return map_list

def TopGraph(file_name, AG, num_q, shortest_length_AG):
    ''' Return the topgraph initial mapping

    Args:
        CX_list (list): the input circuit
        AG (graph): the architecture graph
        num_q (int): the number of qubits in C
    Returns:
        tau (list): the topgraph initial mapping
    '''
    from circuittransform.Li.inimap_c import _tau_bstg_
    #from circuittransform.Li.inimap import _tau_bstg_
    import json
    path = 'C:/ProgramData/Anaconda3/Lib/site-packages/circuittransform/inputs/QASM example/CNOT_list/'
    '''gen CX list'''
    current_path = path + file_name + '.txt'
    with open(current_path, 'r') as f:        
        sqn = json.loads(f.read())
    C = sqn
    '''get map list '''
    tau = _tau_bstg_(C, AG, num_q)
    map_list = [-1] * len(AG.nodes())
    for q_log in range(len(AG.nodes())):
        if q_log in tau.keys():
            map_list[q_log] = tau[q_log]
    '''extend and convert tau to a list where index represents log and value phy'''
# =============================================================================
#     for q in range(num_q):
#         if not q in tau:
#             for index in range(num_q):
#                 if tau[index] == num_q:
#                     tau[index] = q
#                     break
# =============================================================================
    #map_list = MapConvert2(tau, file_name, num_q, shortest_length_AG)
    map_list = FinishMap(map_list)
    return map_list

def WeightedGraph(file_name, AG, num_q, shortest_length_AG):
    ''' Return the WeightedGraph initial mapping

    Args:
        CX_list (list): the input circuit
        AG (graph): the architecture graph
        num_q (int): the number of qubits in C
    Returns:
        tau (list): the topgraph initial mapping
    '''
    #from circuittransform.Li.inimap import _tau_bsg_
    from circuittransform.Li.inimap_c import _tau_bsg_
    import json
    path = 'C:/ProgramData/Anaconda3/Lib/site-packages/circuittransform/inputs/QASM example/CNOT_list/'
    '''gen CX list'''
    current_path = path + file_name + '.txt'
    with open(current_path, 'r') as f:        
        sqn = json.loads(f.read())
    C = sqn
    '''get map list '''
    tau = _tau_bsg_(C, AG)
    map_list = [-1] * len(AG.nodes())
    for q_log in range(len(AG.nodes())):
        if q_log in tau.keys():
            map_list[q_log] = tau[q_log]
    '''extend and convert tau to a list where index represents log and value phy'''
# =============================================================================
#     for q in range(num_q):
#         if not q in tau:
#             for index in range(num_q):
#                 if tau[index] == num_q:
#                     tau[index] = q
#                     break
# =============================================================================
    #map_list = MapConvert2(tau, file_name, num_q, shortest_length_AG)
    map_list = FinishMap(map_list)
    return map_list

if __name__ == '__main__':
    method_AG = ['IBM QX20']
    AG = ct.GenerateArchitectureGraph(20, method_AG)
    '''calculate shortest path and its length'''
    res = ct.ShortestPath(AG)
    shortest_path_AG = res[1]
    shortest_length_AG = (res[0], res[2])
    #map_list = TopGraph('4mod5-v0_20.qasm', AG, 20, shortest_length_AG[0])
    map_list = WeightedGraph('4mod5-v0_20.qasm', AG, 20, shortest_length_AG[0])
    print(map_list)