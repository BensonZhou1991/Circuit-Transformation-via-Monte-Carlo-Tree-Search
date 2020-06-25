# -*- coding: utf-8 -*-
"""
Created on Tue May 14 22:14:39 2019

@author: zxz58
"""
import networkx as nx
from inputs.operationU import OperationCNOT, OperationSingle
from inputs.inputgenerator import GenerateDependency

def QiskitGateToOperation(Gate, flag_single=False):
    '''
    ATTENTION: The input gate should only be cx or single-qubit gates
    convert a Qiskit Gate object to OperationU
    only support CNOT and single qubit gates
    flag_single: whether we should convert single qubit gate
    '''

    '''old qiskit version'''
# =============================================================================
#     if Gate.name == 'cx':
#         qargs = Gate.qargs
#         return OperationCNOT(qargs[0], qargs[1])
# =============================================================================
    '''new qiskit version'''
    if Gate[0].name == 'cx':
        qargs = Gate[1]
        return OperationCNOT(qargs[0], qargs[1])
    else:
        if flag_single == True:
            '''convert single-qbuit gate'''
            qargs = Gate[1]
            if len(qargs) == 1:
                return OperationSingle(q_in=qargs[0],
                                       para=Gate[0].params,
                                       name=Gate[0].name)       
    return None
        

def QiskitCircuitToDG(cir, flag_single=False):
    '''
    convert Qiskit circuit to DG
    support only CNOT and single qubit gates
    flag_single: whether we should convert single qubit gate
    '''
    from operation import OperationToDependencyGraph
    operations = []
    num_unidentified_gates = 0
    qregs = cir.qregs
    if len(qregs) > 1:
        raise Exception('quantum circuit has more than 1 quantum register')
    q = qregs[0]
    data = cir.data
    for gate in data:
        operation = QiskitGateToOperation(gate, flag_single)
        if operation == None:
            num_unidentified_gates += 1
        else:
            operations.append(operation)
    GenerateDependency(operations, q.size)
    DG = OperationToDependencyGraph(operations)
    
    return DG, num_unidentified_gates, q, operations