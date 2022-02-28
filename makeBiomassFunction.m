function makeBiomassFunction
load('D:\Dropbox\Databases\BIGG\iMM904.mat')
model = iND750;
biomassComponents = model.mets(find(model.S(:,find(model.c))));

protein_building_blocks = {'ala__L_c';'arg__L_c';'asn__L_c';'asp__L_c';'cys__L_c';...
    'gln__L_c';'glu__L_c';'gly_c';'his__L_c';'ile__L_c';...
    'leu__L_c';'lys__L_c';'met__L_c';'phe__L_c';'pro__L_c';...
    'ser__L_c';'thr__L_c';'trp__L_c';'tyr__L_c';'val__L_c'};

dna_building_blocks = {'datp_c';'dctp_c';'dgtp_c';'dttp_c'};

rna_building_blocks = {'ctp_c';'gtp_c';'utp_c';'atp_c'};

lipids = {};

diff = setdiff(biomassComponents,[protein_building_blocks;dna_building_blocks;rna_building_blocks])
end