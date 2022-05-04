# -*- coding: utf-8 -*-
"""
Created on Wed May 15 17:31:51 2019

@author: snmen
"""

import sys
import os
import cbmpy

baseFolder = sys.argv[1] #'D:\Dropbox\Databases\BIGG' #sys.argv[1]
inputFileName = sys.argv[2] #'bigg_85_more_fixediNF517_fixed_forMeneco.xml' #'iJO1366.xml' #sys.argv[2]
outputFileName = sys.argv[3] #'bigg_85_more_fixediNF517_fixed_forMeneco_l2.sbml' #sys.argv[3]

os.chdir(baseFolder)
cmod = cbmpy.readSBML3FBC(inputFileName)
cbmpy.CBWrite.writeCOBRASBML(cmod, outputFileName)