# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 11:53:26 2019

@author: Xiangzhen Zhou
"""

import copy
import networkx as nx
import numpy as np
from qiskit import QuantumCircuit, QuantumRegister
import os
from qiskit.extensions import standard
from qiskit.converters import circuit_to_dag
from inputs.operationU import OperationCNOT

def CreateCNOTRandomly(q_reg, num_CNOT, cir = None):
    '''
    generate CNOT operation randomly
    input:
        q_reg: quantum register
        cir: if have, input generated operations to this quantum circuit
    return:
        list of all operatuions
    '''
    
    q = q_reg
    num_qubits = len(q)
    # store all CNOT operation instances
    total_CNOT = []
    # seed for generating random input qubits for each CNOT operations
    seed = np.array(range(num_qubits))
    CNOT_input = np.zeros((num_CNOT, 2))
    for i in range(num_CNOT):
        np.random.shuffle(seed)
        CNOT_input[i, 0:2] = seed[0:2]
    '''generate OperationCNOT instances assuming no swap of CNOTs'''
    # store what is each qubit occupied by currently
    q_occupancy = [[]]*num_qubits
    for i in range(num_CNOT):
        q_c = q[int(CNOT_input[i, 0])]
        q_t = q[int(CNOT_input[i, 1])]
        o_d = []
        if q_occupancy[q_c[1]] != []:
            o_d.append(q_occupancy[q_c[1]])
        if q_occupancy[q_t[1]] != []:
            o_d.append(q_occupancy[q_t[1]])
        new_CNOT = OperationCNOT(q_c, q_t, o_d)
        total_CNOT.append(new_CNOT)
        # refresh q_occupancy
        q_occupancy[q_c[1]] = new_CNOT
        q_occupancy[q_t[1]] = new_CNOT
        
        if isinstance(cir, QuantumCircuit):
            new_CNOT.ConductOperation(cir)
    
    return total_CNOT

def CreateCNOTRandomlyOneLayer(q_log, num_CNOT, AG=None):
    '''
    generate CNOT operation randomly in only one layer
    it may contain executable CNOT gates
    input:
        q_reg: quantum register
        list of all operatuions
    output:
        list [(v_c, v_t), (v_c, v_t), ...]
        list [operation, ...]
    '''
    q = q_log
    num_qubits = len(q)
    pos = list(range(num_qubits))
    np.random.shuffle(pos)
    CNOT_operations = []
    CNOT_list = []
    for i in range(num_CNOT):
        q_c_pos = pos.pop()
        q_t_pos = pos.pop()
        q_c = q[q_c_pos]
        q_t = q[q_t_pos]
        new_CNOT = OperationCNOT(q_c, q_t, [])
        CNOT_operations.append(new_CNOT)
        CNOT_list.append((q_c_pos, q_t_pos))
    return CNOT_list, CNOT_operations

def CreateCircuitFromQASM(file, path=None):
    #QASM_file = open('inputs/QASM example/' + file, 'r')
    if path == None:
        path = './inputs/QASM example/'
    #QASM_file = open('D:/anaconda3/envs/quantum/Lib/site-packages/circuittransform/inputs/QASM example/100_cx_increase_q/' + file,
    #                 'r')
    QASM_file = open(path + file, 'r')
    iter_f = iter(QASM_file)
    QASM = ''
    for line in iter_f: #遍历文件，一行行遍历，读取文本
        QASM = QASM + line
    #print(QASM)
    cir = QuantumCircuit.from_qasm_str(QASM)
    QASM_file.close
    
    return cir

def CreateQASMFilesFromExample():
    path = 'inputs/QASM example/'
    files = os.listdir(path)
    
    return files

def GenerateEdgeofArchitectureGraph(vertex, method):
    edge = []
    num_vertex = len(vertex)
    if method[0] == 'circle' or method[0] == 'directed circle':        
        for i in range(num_vertex-1):
            edge.append((i, i+1))
        edge.append((num_vertex-1, 0))  
    
    '''
    grid architecturew with length = additional_arg[0], width = additional_arg[1]
    '''
    if method[0] == 'grid' or method[0] == 'directed grid':
        length = method[1]
        width = method[2]
        for raw in range(width-1):
            for col in range(length-1):
                current_v = col + raw*length
                edge.append((current_v, current_v + 1))
                edge.append((current_v, current_v + length))
        for raw in range(width-1):
            current_v = (length - 1) + raw*length
            edge.append((current_v, current_v + length))
        for col in range(length-1):
            current_v = col + (width - 1)*length
            edge.append((current_v, current_v + 1))
            
    if method[0] == 'grid2':
        length = method[1]
        width = method[2]
        for raw in range(width-1):
            for col in range(length-1):
                current_v = col + raw*length
                edge.append((current_v, current_v + 1))
                edge.append((current_v, current_v + length))
        for raw in range(width-1):
            current_v = (length - 1) + raw*length
            edge.append((current_v, current_v + length))
        for col in range(length-1):
            current_v = col + (width - 1)*length
            edge.append((current_v, current_v + 1))
        
        current_node1 = length*width
        for raw in [0, width-1]:
            for col in range(length):
                current_node2 = raw*length + col
                edge.append((current_node1, current_node2))
                current_node1 += 1
                
        for raw in [0, length-1]:
            for col in range(width):
                current_node2 = raw + col*length
                edge.append((current_node1, current_node2))
                current_node1 += 1
            
    return edge

def GenerateArchitectureGraph(num_vertex, method, draw_architecture_graph = False):
    '''
    generate architecture graph
    Input:
        method:
            circlr
            grid
            grid2
            IBM QX3
            IBM QX4
            IBM QX5
            IBM QX20
            IBM J-P
            IBM A-B-S
            directed grid
            directed circle
            example in paper
    '''
    if method == ['IBM QX3']:
        G = GenerateArchitectureGraph(16, ['grid', 8, 2])
        G.remove_edges_from([(1, 9),(4, 5)])
        if draw_architecture_graph == True: nx.draw(G, with_labels=True)
        return G

    if method == ['IBM QX5']:
        G = nx.DiGraph()
        vertex = list(range(16))
        edges = [(1,2), (1,0), (2,3), (3,4), (3,14), (5,4), (6,5), (6,11), (6,7),\
                  (7,10), (8,7), (9,8), (9,10), (11,10), (12,5), (12,11), (12,13),\
                  (13,14), (13,4), (15,14), (15,2), (15,0)]
        G.add_nodes_from(vertex)
        G.add_edges_from(edges)
        if draw_architecture_graph == True: nx.draw(G, with_labels=True)
        return G

    if method == ['IBM QX20']:
        G = nx.Graph()
        vertex = list(range(20))
        edges = [(0,1),(1,2),(2,3),(3,4),(0,5),(1,6),(2,7),(3,8),(4,9),(1,7),(2,6),(3,9),(4,8),\
                 (5,6),(6,7),(7,8),(8,9),(5,10),(6,11),(7,12),(8,13),(9,14),(5,11),(6,10),(7,13),(8,12),\
                 (10,11),(11,12),(12,13),(13,14),(10,15),(11,16),(12,17),(13,18),(14,19),(11,17),(12,16),(13,19),(14,18),\
                 (15,16),(16,17),(17,18),(18,19)]
        G.add_nodes_from(vertex)
        G.add_edges_from(edges)
        if draw_architecture_graph == True: nx.draw(G, with_labels=True)
        return G

    if method == ['IBM J-P']:
        G = nx.Graph()
        vertex = list(range(20))
        edges = [(0,1),(1,2),(2,3),(3,4),(0,5),(4,9),
                 (5,6),(6,7),(7,8),(8,9),(5,10),(7,12),(9,14),
                 (10,11),(11,12),(12,13),(13,14),(10,15),(14,19),
                 (15,16),(16,17),(17,18),(18,19)]
        G.add_nodes_from(vertex)
        G.add_edges_from(edges)
        if draw_architecture_graph == True: nx.draw(G, with_labels=True)
        return G

    if method == ['IBM A-B-S']:
        G = nx.Graph()
        vertex = list(range(20))
        edges = [(0,1),(1,2),(2,3),(3,4),(1,6),(3,8),
                 (5,6),(6,7),(7,8),(8,9),(5,10),(7,12),(9,14),
                 (10,11),(11,12),(12,13),(13,14),(11,16),(13,18),
                 (15,16),(16,17),(17,18),(18,19)]
        G.add_nodes_from(vertex)
        G.add_edges_from(edges)
        if draw_architecture_graph == True: nx.draw(G, with_labels=True)
        return G

    if method == ['IBM QX4']:
        G = nx.DiGraph()
        vertex = list(range(5))
        edges = [(1,0), (2,0), (2,1), (2,4), (3,4), (3,2)]
        G.add_nodes_from(vertex)
        G.add_edges_from(edges)
        if draw_architecture_graph == True: nx.draw(G, with_labels=True)
        return G
    
    if method == ['example in paper']:
        G = nx.DiGraph()
        vertex = list(range(6))
        edges = [(1,2), (1,0), (2,3), (3,4), (5,4), (5,2), (5,0)]
        G.add_nodes_from(vertex)
        G.add_edges_from(edges)
        if draw_architecture_graph == True: nx.draw(G, with_labels=True)
        return G        
    
    if method[0] == 'directed grid' or method[0] == 'directed circle':
        G = nx.DiGraph()
    else:
        G = nx.Graph()
    vertex = list(range(num_vertex))
    G.add_nodes_from(vertex)
    edge = GenerateEdgeofArchitectureGraph(vertex, method)
    G.add_edges_from(edge)
    
    if draw_architecture_graph == True: nx.draw(G, with_labels=True)
    
    return G

def GenerateArchitectureGraphError(num_vertex,
                                   method_AG,
                                   method_error,
                                   draw_architecture_graph = False):
    '''
    AG with dicts error_??? to represent error ??? in each node and edge
    '''
    AG = GenerateArchitectureGraph(num_vertex,
                                   method_AG,
                                   draw_architecture_graph)
    error_CX = {}
    error_single = {}
    if method_error[0] == 'test':
        for node in AG.nodes():
            error_single[node] = (np.random.rand()+0.5) * 0.01
        for edge in AG.edges():
            error_CX[node] = (np.random.rand()+0.5) * 0.001
        AG.error_single = error_single
        AG.error_CX = error_CX
    return AG

def CreatePartyMapRandomly(num_qubit, num_CNOT, q_reg):
    '''
    create party map randomly
    '''
    operation_CNOT = CreateCNOTRandomly(q_reg, num_CNOT)
    party_map = np.eye(num_qubit)
    for operation in operation_CNOT:
        c_raw = operation.control_qubit[1]
        t_raw = operation.target_qubit[1]
        party_map[t_raw][:] = np.logical_xor(party_map[t_raw][:], party_map[c_raw][:])
    
    '''set all diagonal elements to 1, it is for testing'''
    '''
    for i in range(num_qubit):
        party_map[i][i] = 1
    '''
    return party_map, operation_CNOT

def GenerateDependency(operations, num_q):
    '''Generate Dependency to operations according to the order'''
    dic = {}
    for i in range(num_q):
        dic[i] = None
    
    for operation in operations:
        qubits = operation.involve_qubits
        for q in qubits:
            if isinstance(q, int):
                q_index = q
            else:
                q_index = q.index
            if dic[q_index] == None:
                dic[q_index] = operation
            else:
                dependent_operation = dic[q_index]
                if not dependent_operation in operation.dependent_operations:
                    operation.dependent_operations.append(dependent_operation)
                dic[q_index] = operation

def CreateDGfromQASMfile(QASM_file, flag_single=False):
    '''
    convert QASM file to cir and DG
    flag_single: whether we should convert single qubit gate
    output:
        circuit, (DG, num_unidentified_gates, quantum register, operations)
    '''
    from Qiskitconverter.QiskitcircuittoDG import QiskitCircuitToDG
    cir = CreateCircuitFromQASM(QASM_file)
    res = QiskitCircuitToDG(cir, flag_single)
    return cir, res

def GenQasmRandomly(num_file, num_qubit, num_cx, path, start_num=0):
    for i in range(num_file):
        full_path = path + str(start_num+i) + '.qasm'
        print('\rGenerating %d of %d file' %(i+1, num_file), end='')
        with open(full_path, 'w') as f:
            f.write('OPENQASM 2.0;\ninclude "qelib1.inc";')
            f.write('\nqreg q[' + str(num_qubit) + '];')
            f.write('\ncreg c[' + str(num_qubit) + '];')
            for _ in range(num_cx):
                cx1 = str(np.random.randint(num_qubit))
                cx2 = str(np.random.randint(num_qubit))
                while cx2 == cx1: cx2 = str(np.random.randint(num_qubit))
                f.write('\ncx q[' + cx1 + '],q[' + cx2 + '];')
        
if __name__ == '__main__':
    start_num = -1
    for num_qubit in [10]:
        for num_cx in [50,100,150,200,250,300,350,400,450,500]:
            for cir_num in range(10):
                start_num +=1
                GenQasmRandomly(num_file=1, num_qubit=num_qubit, num_cx=num_cx,
                                path='temp/', start_num=start_num)
    '''test GenerateArchitectureGraph'''
    '''
    l=4
    w=3
    num_vertex=l*w+2*(l+w)
    G = GenerateArchitectureGraph(num_vertex, ['grid2', l, w], draw_architecture_graph = True)
    '''
    
    '''test GenerateDependency'''
    '''
    num_q = 5
    num_CNOT = 5
    q = QuantumRegister(num_q, 'q')
    cir = QuantumCircuit(q)
    operations = CreateCNOTRandomly(q, num_CNOT, cir)
    GenerateDependency(operations, num_q)
    print(cir.draw())
    DG = ct.OperationToDependencyGraph(operations)
    nx.draw(DG, with_labels=True)
    '''
    
# =============================================================================
#     q = QuantumRegister(4, 'q')
#     cir = QuantumCircuit(q)
#     cir.cx(q[1], q[0])
#     cir.cx(q[1], q[2])
#     cir.cx(q[2], q[3])
#     cir.h(q[3])
#     res = ct.QiskitCircuitToDG(cir)
#     DG = res[0]
#     print(cir.draw())
#     nx.draw(DG, with_labels=True)
# =============================================================================
