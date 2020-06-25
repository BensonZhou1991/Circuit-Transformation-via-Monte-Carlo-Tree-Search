# -*- coding: utf-8 -*-
"""
Created on Mon May 11 01:50:04 2020

@author: zxz58
"""
import numpy as np

num_gates = results['sum total gates 20']
time = results['sum search time 20']

start_num_q = 5
end_num_q = 14
num_cir_q = 10

cir_index = -1
ave_time = []
ave_gate = []
for num_qubit in range(start_num_q, end_num_q+1):
    current_gate = []
    current_time = []
    for cir_num in range(num_cir_q):
        cir_index += 1
        current_gate.append(num_gates[cir_index])
        current_time.append(time[cir_index])
    ave_gate.append(np.average(current_gate))
    ave_time.append(np.average(current_time))
        
        