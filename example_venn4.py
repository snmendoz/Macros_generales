# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 13:00:04 2020

@author: snmen
"""

import venn
#labels = venn.get_labels([
#            range(10),
#            range(5, 15)
#        ], fill=['number'])
#print(labels)
#fig, ax = venn.venn2(labels, names=['list 1', 'list 2'])
#fig.show()


#a = range(10)
#b = range(5,15)
#c = range(3,8)
#d = range(8,7)

#a = ['a', 'b', 'c', 'd']
#b = ['c', 'd', 'e', 'f']
#c = ['e', 'f', 'g', 'h']
#d = ['a', 'e', 'h', 'i']
#
#a_name = 'list 1'
#b_name = 'list 2'
#c_name = 'list 3'
#d_name = 'list 4'

text_file = open("D:\Dropbox\Research_Projects\Atacama\CRG2.txt", "r")
a = text_file.readlines()

text_file2 = open("D:\Dropbox\Research_Projects\Atacama\M. capsulatus.txt", "r")
b = text_file2.readlines()

text_file3 = open("D:\Dropbox\Research_Projects\Atacama\P. putida.txt", "r")
c = text_file3.readlines()

text_file4 = open("D:\Dropbox\Research_Projects\Atacama\B. subtilis.txt", "r")
d = text_file4.readlines()

a_name = 'CRG2'
b_name = 'M. capsulatus'
c_name = 'P. putida'
d_name = 'B. subtilis'

labels = venn.get_labels([a, b, c, d], fill=['number'])
fig, ax = venn.venn4(labels, names=[a_name, b_name, c_name, d_name])
fig.show()