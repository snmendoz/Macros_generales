# -*- coding: utf-8 -*-
"""
Created on Thu Aug 09 13:48:15 2018

@author: snmen
"""

import os
import cbmpy

os.chdir('D:\Dropbox\Databases\BIGG')
#cmod = cbmpy.readSBML3FBC("bigg_85_fixed_forMeneco_cytosol.sbml")
#cbmpy.CBWrite.writeCOBRASBML(cmod, "bigg_85_fixed_forMeneco_cytosol_l2.sbml")

cmod = cbmpy.readSBML3FBC("bigg_85_fixed_forMeneco.sbml")
cbmpy.CBWrite.writeCOBRASBML(cmod, "bigg_85_fixed_forMeneco_l2.sbml")