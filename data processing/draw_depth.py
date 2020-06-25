# -*- coding: utf-8 -*-
"""
Created on Tue May  5 13:55:25 2020

@author: zxz58
"""
import matplotlib.pyplot as plt
import numpy as np

QASM_file = ['rd84_142',
'adr4_197',
'radd_250',
'z4_268',
'sym6_145',
#'misex1_241',
'rd73_252',
'cycle10_2_110',
'square_root_7',
'sqn_258',
'rd84_253']


width = 2
length = 2

gridspec_kw = {'hspace':0.05}
fig, axs = plt.subplots(length, width, gridspec_kw=gridspec_kw)  # a figure with a 2x2 grid of Axes

for i in range(length):
    for j in range(width):
        name = QASM_file[i*width + j]
        best_index = np.argmin(results[name]['total gates number'])
        depth = results[name]['depth'][best_index]
        dep_min = depth[0]
        dep_ave = depth[1]
        dep_max = depth[2]
        axs[i][j].plot(dep_min, label='Min.')
        axs[i][j].plot(dep_ave, label='Ave.')
        axs[i][j].plot(dep_max, label='Max.')
        axs[i][j].set_title(name)
        axs[i][j].set_xlabel('Decision')
        axs[i][j].set_ylabel('Search Depth')
        axs[i][j].legend(loc='upper right')
fig.savefig('depth.pdf' ,format='pdf', papertype='a4')
            
# =============================================================================
#     search_tree.DrawLogData()
#     plt.plot(search_depth_ave, label='average')
#     plt.plot(search_depth_max, label='maximum')
#     plt.plot(search_depth_min, label='minimum')
#     plt.legend()
#     plt.grid()
#     plt.title('search depth')
#     plt.show()
# =============================================================================
