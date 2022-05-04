# -*- coding: utf-8 -*-
"""
#WORKS WITH PYTHON 2

Created on Thu Jul 18 20:20:31 2019

@author: snmen
"""

import cobra
import os
import sys
#inputName = sys.argv[1]
#outputName = sys.argv[2]
#folder = sys.argv[3]


inputName = 'st_v5_withBounds.xml'
outputName ='st_v5_withBounds.json'
folder = 'D:/Dropbox/Research_Projects/Lactobacillus_bulgaricus/L_bulgaricus/scripts'

#inputName = 'S_thermophillus_v1_03.xml'
#outputName ='S_thermophillus_v1_03.json'
#folder = 'D:\Dropbox\Research_Projects\Streptococcus thermophilus'

#inputName = 'lbu_casein.xml'
#outputName = 'lbu_casein.json'
#folder = 'D:/Dropbox/Research_Projects/Lactobacillus_bulgaricus/L_bulgaricus/scripts'

#inputName = 'LBUL_v1_07_b_lac.xml'
#outputName ='LBUL_v1_07_b_lac.json'
#folder = 'D:/Dropbox/Research_Projects/Lactobacillus_bulgaricus/L_bulgaricus/reconstructions/Manual'

#inputName = 'lb_lactate_weird.xml'
#outputName ='lb_lactate_weird.json'
#folder = 'D:/Dropbox/Research_Projects/Lactobacillus_bulgaricus/L_bulgaricus/scripts'


#inputName = 'plantarumMetadraft.xml'
#outputName ='plantarumMetadraft.json'
#folder = 'D:\Dropbox\Research_Projects\crispatus'


#inputName = 'plantarumMetadraft_bounds.xml'
#outputName ='plantarumMetadraft_bounds.json'
#folder = 'D:\Dropbox\Research_Projects\crispatus'

#inputName = 'yeast841_biomass_pizarro_2007_30_degrees.xml'
#outputName ='yeast841_biomass_pizarro_2007_30_degrees.json'
#folder = 'D:\Dropbox\Research_Projects\protocol_wine_fba'

#inputName = 'LBUL_v1_03_polypeptides3.xml'
#outputName ='LBUL_v1_03_polypeptides3.json'
#folder = 'D:/Dropbox/Research_Projects/Lactobacillus_bulgaricus/L_bulgaricus/reconstructions/Manual'

#inputName = 'LBUL_v1_03_polypeptides4_2.xml'
#outputName ='LBUL_v1_03_polypeptides4_2.json'
#folder = 'D:/Dropbox/Research_Projects/Lactobacillus_bulgaricus/L_bulgaricus/reconstructions/Manual'

#inputName = 'LBUL_v1_03_polypeptides4_1.xml'
#outputName ='LBUL_v1_03_polypeptides4_1.json'
#folder = 'D:/Dropbox/Research_Projects/Lactobacillus_bulgaricus/L_bulgaricus/reconstructions/Manual'

os.chdir(folder)
model = cobra.io.read_sbml_model(inputName)
print(model)
cobra.io.save_json_model(model, outputName)