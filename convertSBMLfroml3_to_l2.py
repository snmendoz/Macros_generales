# -*- coding: utf-8 -*-
"""
Created on Thu Aug 09 10:54:38 2018

@author: snmen
"""
import os
import cbmpy

baseFolder = 'D:\Dropbox\Databases\AuremeInputs' 
os.chdir(baseFolder)

modelsList = open("currentModels.txt").read().splitlines()

for modelName in modelsList:
    os.chdir(os.path.join(baseFolder,modelName))
    cmod = cbmpy.readSBML3FBC("metabolic_model_l3.sbml")
    cbmpy.CBWrite.writeCOBRASBML(cmod, "metabolic_model_cobra.xml")
    os.chdir('D:\Dropbox\Databases\AuremeInputs')
