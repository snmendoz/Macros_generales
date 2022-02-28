function [genes, sequences, geneNames, locus_tag] = readFastaProteinsGeneral(nombreArchivo)

fileID = fopen(nombreArchivo,'r+');
contador=0;
genes=cell(10000,1);
sequences=cell(10000,1);
geneNames=cell(10000,1);
locus_tag=cell(10000,1);

tline=fgetl(fileID);
while ischar(tline)
    
    if strcmp(tline(1), '>')
        info = regexp(tline,' ','split');
        contador = contador + 1; 
        gene = info{1}(2:end);
        genes{contador} = gene;
        if ~isempty(find(cellfun(@isempty,strfind(info,'locus_tag'))==0))
            locus_tag{contador} = regexprep(info{find(cellfun(@isempty,strfind(info,'locus_tag'))==0)},{'[locus_tag=',']'},{'',''});
        end
        tline=fgetl(fileID);
        secuencia = '';
        while ischar(tline) && isempty(strfind(tline,'>'))
            secuencia = [secuencia tline]; 
            tline=fgetl(fileID);
        end
        
        sequences{contador} = secuencia;
    end
    
end

genes = genes(1:contador);
sequences = sequences(1:contador);

b1=find(cellfun(@isempty, genes));
b2=find(cellfun(@isempty, sequences));

if ~isempty(b1) || ~isempty(b2)
    warning('genes or sequences missing')
end
fclose(fileID);

end