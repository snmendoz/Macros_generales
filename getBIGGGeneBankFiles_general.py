# -*- coding: utf-8 -*-
"""
Created on Wed Mar 07 14:43:41 2018

@author: snmen
"""
import os
from Bio import Entrez
import requests

r = requests.get('http://bigg.ucsd.edu/api/v2/models')
r.json()
id = r.json()["results"]
names = []

for i in id:
    name = i["bigg_id"]
    names.append(name)

os.chdir('D:/Dropbox/Databases/BIGG/genomes/all')
nested_list = []
l = []
l.append('ModelID')
l.append('GenomeID')
nested_list.append(l)

for model in names:
    r = requests.get('http://bigg.ucsd.edu/api/v2/models/' + model)
    r.json()
    id = r.json()["genome_name"]
    
    if id:

        Entrez.email = "s.n.mendozafarias@vu.nl"
        handle = Entrez.esearch(db="nuccore", term=id)
        record = Entrez.read(handle)
        gi_list = record["IdList"]
        print(gi_list)
    
        gi_str = ",".join(gi_list)
        handle = Entrez.efetch(db="nuccore", id=gi_str, rettype="fasta_cds_aa", retmode="text")
    
        text = handle.read()
        #print(text)
        file = open(id + "_full.faa","w")
        file.write(text) 
        file.close() 
        if os.path.isfile(model + "_full.faa"): 
            os.remove(model + "_full.faa") 
            os.rename(id + "_full.faa", model + "_full.faa")
        else:
            os.rename(id + "_full.faa", model + "_full.faa")
        
        handle = Entrez.efetch(db="nuccore", id=gi_str, rettype="gbwithparts", retmode="text")
        text = handle.read()
        #print(text)
        file = open(id + ".gb","w")
        file.write(text) 
        file.close() 
        if os.path.isfile(model + ".gb"): 
            os.remove(model + ".gb") 
            os.rename(id + ".gb", model + ".gb")
        else:
            os.rename(id + ".gb", model + ".gb")
        print(model)
        print(id)
        
        l = []
        l.append(model)
        l.append(id)
        nested_list.append(l)

    else:
        l = []
        l.append(model)
        l.append('null')
        nested_list.append(l)
            
with open('genomeIDs.txt', 'w') as file:
    file.writelines('\t'.join(i) + '\n' for i in nested_list)

#handle = Entrez.efetch(db="nuccore", id="NC_009004.1", rettype="gb", retmode="text")
#text = handle.read()
#print(text)

#records = SeqIO.parse(handle, "gb")
#for record in records:
#    print("%s, length %i, with %i features" % (record.name, len(record), len(record.features)))
    