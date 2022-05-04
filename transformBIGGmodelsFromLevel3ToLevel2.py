# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 16:30:51 2018

@author: snmen
"""
import os
import cbmpy

os.chdir('D:\Dropbox\Databases\BIGG')

with open('modelNames.txt') as f:
    modelNames = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each line
modelNames = [x.strip() for x in modelNames] 

excluded = {'iAT_PLT_636'}

for modelName in content:
    if modelName not in excluded:
        cmod = cbmpy.readSBML3FBC(modelName + '.xml')
        os.chdir('D:\Dropbox\Databases\BIGG\SBML_L2')
        cbmpy.CBWrite.writeCOBRASBML(cmod, modelName + '_l2.xml')
        os.chdir('D:\Dropbox\Databases\BIGG')
    
