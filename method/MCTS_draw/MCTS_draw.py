# -*- coding: utf-8 -*-
"""
Created on Fri May  8 00:08:40 2020

@author: zxz58
"""
from method.MCtree import MCTree
from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister
from operation import FindExecutableNode, ExecuteAllPossibileNodesInDG
from inputs.inputgenerator import CreateDGfromQASMfile

def AddSwap(cir, swap, qreg):
    cir.cx(qreg[swap[0]], qreg[swap[1]])
    cir.cx(qreg[swap[1]], qreg[swap[0]])
    cir.cx(qreg[swap[0]], qreg[swap[1]])

def MapToNaive(mapping, cir, qreg):
    current_list = mapping.MapToList()
    num_q_log = len(current_list)
    for i in range(num_q_log):
        if current_list[i] != i:
            for j in range(i+1, num_q_log):
                if current_list[j] == i:
                    swap = (current_list[i], current_list[j])
                    mapping.RenewMapViaExchangeCod(swap[0], swap[1])
                    AddSwap(cir, swap, qreg)
                    current_list = mapping.MapToList()
                    break
            if current_list[i] != i:
                swap = (current_list[i], i)
                mapping.RenewMapViaExchangeCod(swap[0], swap[1])
                AddSwap(cir, swap, qreg)
                current_list = mapping.MapToList()
    print(mapping.MapToList())

def AddOperationInCircuit(op, cir, qreg, mapping):
    qubits_log = op.involve_qubits_list
    qubits = []
    for q_log in qubits_log:
        q_phy = mapping.LogToPhy(q_log)
        qubits.append(q_phy)
    gate_name = op.name
    if gate_name == 'CX':
        cir.cx(qreg[qubits[0]], qreg[qubits[1]])
        return None
    if gate_name == 't':
        cir.t(qreg[qubits[0]])
        return None
    if gate_name == 'h':
        cir.h(qreg[qubits[0]])
        return None
    if gate_name == 'tdg':
        cir.tdg(qreg[qubits[0]])
        return None
    if gate_name == 'x':
        cir.x(qreg[qubits[0]])
        return None
    if gate_name == 'rx':
        cir.rx(float(op.para[0]), qreg[qubits[0]])
        return None
    if gate_name == 'ry':
        cir.ry(float(op.para[0]), qreg[qubits[0]])
        return None
    if gate_name == 'rz':
        cir.rz(float(op.para[0]), qreg[qubits[0]])
        return None
    raise(Exception('Unkonow gate', gate_name))

def MCTSToQiskitCir(MCTree, file, convert_final_map=False):
    '''
    create physical circuit with all single qubits with a MCTree
    '''
    AG = MCTree.AG
    res_qasm = CreateDGfromQASMfile(file, flag_single=True)
    DG = res_qasm[1][0]
    '''initialize physical quantum circuit'''
    q_phy = QuantumRegister(len(AG.nodes()), 'q')
    c_phy = ClassicalRegister(len(AG.nodes()), 'c')
    cir_phy = QuantumCircuit(q_phy)
    cir_phy.add_register(c_phy)
    '''cope with node 0'''
    current_node = 0
    current_map = MCTree.nodes[current_node]['mapping']
    executed_vertex = []
    executable_vertex = FindExecutableNode(DG)
    res = ExecuteAllPossibileNodesInDG(executable_vertex, executed_vertex, AG, 
                                       DG, current_map,
                                       out_removed_nodes=True,
                                       draw=False, DiG=None, edges_DiG=None)
    executed_vertex, executable_vertex, executed_vertex_current = res
    for vertex in executed_vertex_current:
        op = DG.nodes[vertex]['operation']
        AddOperationInCircuit(op, cir_phy, q_phy, current_map)     
    while MCTree.nodes[current_node]['num_remain_vertex'] != 0:
        sons = list(MCTree.successors(current_node))
        current_node = sons[0]
        swap = MCTree.nodes[current_node]['added_SWAP'][0]
        AddSwap(cir_phy, swap, q_phy)
        current_map = MCTree.nodes[current_node]['mapping']
        res = ExecuteAllPossibileNodesInDG(executable_vertex, executed_vertex, AG, 
                                           DG, current_map,
                                           out_removed_nodes=True,
                                           draw=False, DiG=None, edges_DiG=None)
        executed_vertex, executable_vertex, executed_vertex_current = res
        for vertex in executed_vertex_current:
            op = DG.nodes[vertex]['operation']
            AddOperationInCircuit(op, cir_phy, q_phy, current_map)
    # convert the final map to naive one without considering AG
    if convert_final_map == True:
        MapToNaive(current_map, cir_phy, q_phy)
    return cir_phy


MCTree.ToQiskitCir = MCTSToQiskitCir