# -*- coding: utf-8 -*-
"""
Created on Wed May 22 10:31:16 2019

@author: snmen
"""

import venn
import sys
import os
import matplotlib.pyplot as plt

baseFolder = sys.argv[1]

os.chdir(baseFolder)

groupNames = open("groupNames.txt").read().splitlines()
list1 = open("list1.txt").read().splitlines()
list2 = open("list2.txt").read().splitlines()
list3 = open("list3.txt").read().splitlines()
list4 = open("list4.txt").read().splitlines()
list5 = open("list5.txt").read().splitlines()

labels = venn.get_labels([list1,list2,list3,list4,list5], fill=['number', 'str'])
fig, ax = venn.venn5(labels, names = groupNames)
fig.savefig("comparison.pdf", bbox_inches='tight')

