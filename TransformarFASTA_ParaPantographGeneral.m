function [Genes, Secuencias, oldGenes, Enzimas, EC]=TransformarFASTA_ParaPantographGeneral(nombreArchivo, nombreNuevoArchivo, iniciador, writeFile, considerFinalGene)

if nargin <3 || isempty(iniciador)
    iniciador = 1;
end
if nargin <4 || isempty(writeFile)
    writeFile = 1;
end

if nargin <5 || isempty(considerFinalGene)
    considerFinalGene = 0;
end
    
if iniciador == 1
    principioGen = '    gene           ';
elseif iniciador ==2
    principioGen = '    CDS            ';
end

if considerFinalGene
    finalGene = '    CDS            ';
end

fileID = fopen(nombreArchivo,'r+');
Genes=cell(10000,1);
oldGenes=cell(10000,1);
Enzimas=cell(10000,1);
EC=cell(10000,1);
Secuencias=cell(10000,1);

contador=0;
locus_tag = 'None';
old_locus_tag = 'None';
Translation = 'None';
enzyme = 'None';
ECnumber = 'None';

tline=fgetl(fileID);
while ischar(tline)
%     if ~isempty(findstr(tline, 'ORIGIN'))
%         break;
%     end
    
    if ~isempty(findstr(tline, principioGen))
        if contador >=1 && length(find(cellfun(@isempty, Genes(1:contador))))==1
            contador = contador - 1;
        end
        contador = contador + 1;
        locus_tag = '';
        old_locus_tag = '';
        Translation = '';
        enzyme = '';
        ECnumber = '';
    end
    
    if considerFinalGene && ~isempty(findstr(tline, finalGene))
        if isempty(locus_tag); locus_tag = 'None'; end
        if isempty(old_locus_tag); old_locus_tag = 'None'; end
        if isempty(Translation); Translation = 'None'; end
        if isempty(enzyme); enzyme = 'None'; end 
        if isempty(ECnumber); ECnumber = 'None'; end
    end
    
    if ~isempty(findstr(tline, '/locus_tag="')) && isempty(locus_tag) %|| ~isempty(findstr(tline, '/gene="')))
        split=regexp(tline,'\"','split');
        locus_tag=split(end-1);
        Genes(contador) = locus_tag;
    end
    
    if ~isempty(findstr(tline, '/old_locus_tag="')) && isempty(old_locus_tag)
        split=regexp(tline,'\"','split');
        old_locus_tag=split(end-1);
        oldGenes(contador) = old_locus_tag;
    end
    
    if ~isempty(findstr(tline, '/translation="')) && isempty(Translation)
        split=regexp(tline,'\"','split');
        if length(strfind(tline,'"'))==2
            
            Translation=split{end-1};
        else
            Translation=split{end};
            tline=fgetl(fileID);
            while isempty(findstr(tline, '"'))
                Translation=[Translation regexprep(tline,' ','')];
                Translation=regexprep(Translation,'\t','');
                Translation=regexprep(Translation,'\n','');
                tline=fgetl(fileID);
            end
            split=regexp(tline,'\"','split');
            Translation=[Translation regexprep(split{1},' ','')];
        end
        Translation=regexprep(Translation,'\t','');
        Translation=regexprep(Translation,'\n','');
        Secuencias(contador)={Translation};
    end
    
    if (~isempty(findstr(tline, '/product="')) || ~isempty(findstr(tline, '/function="'))) && isempty(enzyme)
        if length(strfind(tline,'"'))==2
            split=regexp(tline,'\"','split');
            enzyme = split(end-1);
            Enzimas(contador)=enzyme;
        else
            split=regexp(tline,'\"','split');
            Enzimas(contador)=split(end);
            tline=fgetl(fileID);
            split=regexp(tline,'\"','split');
            Enzimas{contador}=[Enzimas{contador} ' ' regexprep(split{1},' ','')];
            enzyme = Enzimas(contador);
        end
    end
    
    if (~isempty(findstr(tline, '/EC_number="')) || ~isempty(findstr(tline, '/EC="'))) && isempty(ECnumber)
        split=regexp(tline,'\"','split');
        ECnumber = split(end-1);
        EC(contador)= ECnumber;
    end
    
    tline=fgetl(fileID);
end
fclose(fileID);
% xlswrite(['Proteinas' nombreArchivo],[Genes, Enzimas, EC ])

Genes=Genes(1:contador);
oldGenes=oldGenes(1:contador);
EC=EC(1:contador);
Enzimas=Enzimas(1:contador);
Secuencias=Secuencias(1:contador);
% Enzimas=Enzimas(1:pos-1);
% EC=EC(1:pos-1);
if writeFile
    if ~isempty(strfind(nombreNuevoArchivo,'.'))
        fileID2 = fopen(nombreNuevoArchivo,'wt');
    else
        fileID2 = fopen([nombreNuevoArchivo '.fasta'],'wt');
    end
    for i=1:length(Genes)
        if ~isempty(Secuencias{i})
            fprintf(fileID2, '%s\n', ['>' Genes{i}])
            fprintf(fileID2, '%s\n', Secuencias{i})
        end
    end
    fclose(fileID2);
end
end