# -*- coding: utf-8 -*-
"""
Created on Wed May 13 17:57:27 2020

@author: zxz58
"""

'''
this module is to get the number of CNOT gates in QASM files
'''

# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 23:21:22 2020

@author: zxz58
"""
import circuittransform as ct
   
bench_mark = 'Cambridge' #'Q20', 'Cambridge', 'Benchmark'

'''Q20 Benchmark1'''
if bench_mark == 'Q20':
    QASM_file = ['4mod5-v1_22.qasm',
    'mod5mils_65.qasm',
    'alu-v0_27.qasm',
    'decod24-v2_43.qasm',
    '4gt13_92.qasm',
    'ising_model_10.qasm',
    'ising_model_13.qasm',
    'ising_model_16.qasm',
    'qft_10.qasm',
    'qft_16.qasm',
    'rd84_142.qasm',
    'adr4_197.qasm',
    'radd_250.qasm',
    'z4_268.qasm',
    'sym6_145.qasm',
    'misex1_241.qasm',
    'rd73_252.qasm',
    'cycle10_2_110.qasm',
    'square_root_7.qasm',
    'sqn_258.qasm',
    'rd84_253.qasm',
    'co14_215.qasm',
    'sym9_193.qasm']

'''MCTS papaer benchmark'''
if bench_mark == 'Benchmark':
    QASM_file = ['rd84_142.qasm',
    'adr4_197.qasm',
    'radd_250.qasm',
    'z4_268.qasm',
    'sym6_145.qasm',
    'misex1_241.qasm',
    'rd73_252.qasm',
    'cycle10_2_110.qasm',
    'square_root_7.qasm',
    'sqn_258.qasm',
    'rd84_253.qasm']

'''Q20 Benchmark2'''
if bench_mark == 'Cambridge':
    QASM_file = ['graycode6_47.qasm',
    'xor5_254.qasm',
    'ex1_226.qasm',
    '4gt11_84.qasm',
    'ex-1_166.qasm',
    'ham3_102.qasm',
    '4mod5-v0_20.qasm',
    '4mod5-v1_22.qasm',
    'mod5d1_63.qasm',
    '4gt11_83.qasm',
    '4gt11_82.qasm',
    'rd32-v0_66.qasm',
    'mod5mils_65.qasm',
    '4mod5-v0_19.qasm',
    'rd32-v1_68.qasm',
    'alu-v0_27.qasm',
    '3_17_13.qasm',
    '4mod5-v1_24.qasm',
    'alu-v1_29.qasm',
    'alu-v1_28.qasm',
    'alu-v3_35.qasm',
    'alu-v2_33.qasm',
    'alu-v4_37.qasm',
    'miller_11.qasm',
    'decod24-v0_38.qasm',
    'alu-v3_34.qasm',
    'decod24-v2_43.qasm',
    'mod5d2_64.qasm',
    '4gt13_92.qasm',
    '4gt13-v1_93.qasm',
    'one-two-three-v2_100.qasm',
    '4mod5-v1_23.qasm',
    '4mod5-v0_18.qasm',
    'one-two-three-v3_101.qasm',
    '4mod5-bdd_287.qasm',
    'decod24-bdd_294.qasm',
    '4gt5_75.qasm',
    'alu-v0_26.qasm',
    'rd32_270.qasm',
    'alu-bdd_288.qasm',
    'decod24-v1_41.qasm',
    '4gt5_76.qasm',
    '4gt13_91.qasm',
    '4gt13_90.qasm',
    'alu-v4_36.qasm',
    '4gt5_77.qasm',
    'one-two-three-v1_99.qasm',
    'rd53_138.qasm',
    'one-two-three-v0_98.qasm',
    '4gt10-v1_81.qasm',
    'decod24-v3_45.qasm',
    'aj-e11_165.qasm',
    '4mod7-v0_94.qasm',
    'alu-v2_32.qasm',
    '4mod7-v1_96.qasm',
    'cnt3-5_179.qasm',
    'mod10_176.qasm',
    '4gt4-v0_80.qasm',
    '4gt12-v0_88.qasm',
    '0410184_169.qasm',
    '4_49_16.qasm',
    '4gt12-v1_89.qasm',
    '4gt4-v0_79.qasm',
    'hwb4_49.qasm',
    '4gt4-v0_78.qasm',
    'mod10_171.qasm',
    '4gt12-v0_87.qasm',
    '4gt12-v0_86.qasm',
    '4gt4-v0_72.qasm',
    '4gt4-v1_74.qasm',
    'mini-alu_167.qasm',
    'one-two-three-v0_97.qasm',
    'rd53_135.qasm',
    'ham7_104.qasm',
    'decod24-enable_126.qasm',
    'mod8-10_178.qasm',
    '4gt4-v0_73.qasm',
    'ex3_229.qasm',
    'mod8-10_177.qasm',
    'alu-v2_31.qasm',
    'C17_204.qasm',
    'rd53_131.qasm',
    'alu-v2_30.qasm',
    'mod5adder_127.qasm',
    'rd53_133.qasm',
    'majority_239.qasm',
    'ex2_227.qasm',
    'cm82a_208.qasm',
    'sf_276.qasm',
    'sf_274.qasm',
    'con1_216.qasm',
    'rd53_130.qasm',
    'f2_232.qasm',
    'rd53_251.qasm',
    'hwb5_53.qasm',
    'radd_250.qasm',
    'rd73_252.qasm',
    'cycle10_2_110.qasm',
    'hwb6_56.qasm',
    'cm85a_209.qasm',
    'rd84_253.qasm',
    'root_255.qasm',
    'mlp4_245.qasm',
    'urf2_277.qasm',
    'sym9_148.qasm',
    'hwb7_59.qasm',
    'clip_206.qasm',
    'sym9_193.qasm',
    'dist_223.qasm',
    'sao2_257.qasm',
    'urf5_280.qasm',
    'urf1_278.qasm',
    'sym10_262.qasm',
    'hwb8_113.qasm',
    #'test3.qasm'
    ]
cx_count = []

for i in range(len(QASM_file)):
    file = QASM_file[i]
    if file[-5:] != '.qasm': file += '.qasm'
    res_qasm = ct.CreateDGfromQASMfile(file)
    name = file[0:-5]
 
    '''generate dependency graph'''
    DG = res_qasm[1][0]
    cx_count.append(len(DG.nodes()))