# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 17:49:45 2019

@author: snmen
"""

# Import the library
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
 
# Custom text labels: change the label of group A
#                 (1 , 2  ,1-2, 3  ,1-3,2-3,1-2-3)
v=venn3(subsets = (90, 380, 85, 779,518,74 ,482), set_labels = ('Microbacterium CRG1', 'M. tuberculosis', 'S. coelicolor'))
#v.get_label_by_id('A').set_text('Microbacterium')
plt.show()

#                 (1 , 2  ,1-2, 3  ,1-3,2-3,1-2-3)
v=venn3(subsets = (98, 380, 85, 778,519,78,478), set_labels = ('Microbacterium CRG2', 'M. tuberculosis', 'S. coelicolor'))
#v.get_label_by_id('A').set_text('Microbacterium')
plt.show()
 
## Line style: can be 'dashed' or 'dotted' for example
#v=venn3(subsets = (10, 8, 22, 6,9,4,2), set_labels = ('Group A', 'Group B', 'Group C'))
#c=venn3_circles(subsets = (10, 8, 22, 6,9,4,2), linestyle='dashed', linewidth=1, color="grey")
#plt.show()
# 
## Change one group only
#v=venn3(subsets = (10, 8, 22, 6,9,4,2), set_labels = ('Group A', 'Group B', 'Group C'))
#c=venn3_circles(subsets = (10, 8, 22, 6,9,4,2), linestyle='dashed', linewidth=1, color="grey")
#c[0].set_lw(8.0)
#c[0].set_ls('dotted')
#c[0].set_color('skyblue')
#plt.show()
# 
## Color
#v.get_patch_by_id('100').set_alpha(1.0)
#v.get_patch_by_id('100').set_color('white')
#plt.show()
