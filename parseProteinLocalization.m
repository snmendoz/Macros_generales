function sublocalization = parseProteinLocalization(fileName,type)

switch type
    case 'deeploc'
        info2 = readCSVFile(fileName, ',');
        sublocalization = info2(2:end,1:2);
    case 'wolf'
        fid = fopen(fileName,'r+');
        tline = fgetl(fid);
%         info = regexp(tline, '(XP_[^ ]* details [^ ]*:)','tokens')
        seqs = regexp(tline, '(DEHA2[^ ]*) ','tokens');
        loc = regexp(tline, 'DEHA2[^ ]* details ([^ ]*):','tokens');
        sublocalization = [cellfun(@(x) x, seqs)' cellfun(@(x) x, loc)'];
        sublocalization = regexprep(sublocalization,{'cysk','cyto','E.R','extr','golg','mito','nucl','plas','pero','vacu','chlo'},{'cytoskeleton','cytosol','endoplasmic_reticulum','extracellular','golgi','microchondrion','nucleus','plasma_membrane','peroxisome','vacuolar_membrane','thylakoid_lumen'});
    case 'busca'
        info2 = readCSVFile_general(fileName, ',');
        sublocalization = [info2(2:end,1) regexprep(info2(2:end,3),'C:','')];
    case 'cello'
        info2 = readCSVFile_general(fileName, '\t');
        sublocalization = [info2(2:end,21) info2(2:end,20)];
        
        
        
end


end