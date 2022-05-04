# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 12:52:48 2019

@author: snmen
"""

import sys
from bioservices.kegg import KEGG
from collections import Counter
import csv
import operator
import os

#fileName = sys.argv[1]
fileName = 'D:\Dropbox\Research_Projects\Review_reconstruction\comparison\ppu\metList.txt'
userNames  = metList = open(fileName).read().splitlines()
#userNames = ['D-Fructose']
#userNames = ['Glycine',
#'(S)-Malate',
#'Iron',
#'L-Tryptophan',
#'Uracil',
#'Urea',
#'Oxygen',
#'D-Fructose',
#'gamma-butyrobetaine',
#'D-Glutamate',
#'D-Glucarate',
#'D-Glucose']
#userNames = ['Pseudomonas Common O Polysaccharide']

k = KEGG()
#data = k.get('C00299')
#print(data)
allIds = []
allUserNames = []
for userName in userNames:
    print(userName)
    data2 = k.find('compound', userName)
    if type(data2) is unicode:
        if data2[0]!='\n':
            dictionary = data2.split('\n'); 
            dictionary.pop()
            possibleIDs = []
            
            wasFound = 0
            for pairs in dictionary:
                #print(userName)
#                print(pairs)
                info = pairs.split('\t')
                id = info[0]
                names = info[1].split('; ')
                pos =  [i for i,x in enumerate(names) if x == userName]
                if pos:
                    possibleIDs.append(id.replace('cpd:',''))
                    wasFound = 1
            if wasFound:
                ids = ';'.join(possibleIDs)
                allIds.append(ids)
                allUserNames.append(userName)
    
finalInfo = dict(zip(allUserNames, allIds))    
#print(userNames)
#print(allIds)
#print(finalInfo)
    
with open('idMapping.csv','wb') as out:
    csv_out=csv.writer(out, delimiter='\t')
    for key in finalInfo.keys():
        csv_out.writerow([key, finalInfo[key]])    