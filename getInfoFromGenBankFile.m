function info = getInfoFromGenBankFile(fileName, starter,considerFinalGene)

if nargin <2 || isempty(starter)
    starter = 1;
end

if nargin <3 || isempty(considerFinalGene)
    considerFinalGene = 0;
end
    
if starter == 1
    principioGen = '    gene           ';
elseif starter ==2
    principioGen = '    CDS            '; 
elseif starter ==3
    principioGen = '     CDS             '; %Christian Hansen format
end

if considerFinalGene
    finalGene = '    CDS            ';
end

fileID = fopen(fileName,'r+');
genes=cell(10000,1);
oldGenes=cell(10000,1);
enzymes=cell(10000,1);
products = cell(10000,1);
functions =cell(10000,1);
EC=cell(10000,1);
sequences=cell(10000,1);
pseudoGenes=cell(10000,1);
protein_ids = cell(10000,1);
names = cell(10000,1);


counter=0;
locus_tag = 'None';
old_locus_tag = 'None';
translation = 'None';
enzyme = 'None';
product = 'None';
function_= 'None';
ECnumber = 'None';
pseudo = '0';
protein_id = 'None';
name = 'None';

tline=fgetl(fileID);
while ischar(tline)
%     if ~isempty(findstr(tline, 'ORIGIN'))
%         break;
%     end
    
    if ~isempty(findstr(tline, principioGen))
        if counter >=1 && length(find(cellfun(@isempty, genes(1:counter))))==1
            counter = counter - 1;
        end
        counter = counter + 1;
        locus_tag = '';
        old_locus_tag = '';
        translation = '';
        enzyme = '';
        ECnumber = '';
        pseudo = '';
        product = '';
        function_= '';
        protein_id = '';
        name = '';
    end
    
    if considerFinalGene && ~isempty(findstr(tline, finalGene))
        if isempty(locus_tag); locus_tag = 'None'; end
        if isempty(old_locus_tag); old_locus_tag = 'None'; end
        if isempty(translation); translation = 'None'; end
        if isempty(enzyme); enzyme = 'None'; end 
        if isempty(ECnumber); ECnumber = 'None'; end
        if isempty(pseudo); pseudo = '0'; end
        if isempty(product); product = 'None'; end
        if isempty(function_); function_ = 'None'; end
        if isempty(protein_id); protein_id = 'None'; end
        if isempty(name); name = 'None'; end

    end
    
    if ~isempty(findstr(tline, '/locus_tag="')) && isempty(locus_tag) %|| ~isempty(findstr(tline, '/gene="')))
        split=regexp(tline,'\"','split');
        locus_tag=split(end-1);
        genes(counter) = locus_tag;
    end
    
    if ~isempty(findstr(tline, '/Name="')) && isempty(name)
        split=regexp(tline,'\"','split');
        name=split(end-1);
        names(counter) = name;
    end
    
    if ~isempty(findstr(tline, '/old_locus_tag="')) && isempty(old_locus_tag)
        split=regexp(tline,'\"','split');
        old_locus_tag=split(end-1);
        oldGenes(counter) = old_locus_tag;
    end
    
    if ~isempty(findstr(tline, '/protein_id="')) && isempty(protein_id)
        split=regexp(tline,'\"','split');
        protein_id=split(end-1);
        protein_ids(counter) = protein_id;
    end
    
    if ~isempty(findstr(tline, '/translation="')) && isempty(translation)
        split=regexp(tline,'\"','split');
        if length(strfind(tline,'"'))==2
            
            translation=split{end-1};
        else
            translation=split{end};
            tline=fgetl(fileID);
            while isempty(findstr(tline, '"'))
                translation=[translation regexprep(tline,' ','')];
                translation=regexprep(translation,'\t','');
                translation=regexprep(translation,'\n','');
                tline=fgetl(fileID);
            end
            split=regexp(tline,'\"','split');
            translation=[translation regexprep(split{1},' ','')];
        end
        translation=regexprep(translation,'\t','');
        translation=regexprep(translation,'\n','');
        sequences(counter)={translation};
    end
    
    if ~isempty(findstr(tline, '/product="')) && isempty(product)
        if length(strfind(tline,'"'))==2
            split=regexp(tline,'\"','split');
            product = split(end-1);
            products(counter)=product;
        else
            split=regexp(tline,'\"','split');
            products(counter)=split(end);
            tline=fgetl(fileID);
            split=regexp(tline,'\"','split');
            products{counter}=[products{counter} ' ' regexprep(split{1},' ','')];
            product = enzymes(counter);
        end
    end
    
    if ~isempty(findstr(tline, '/function="')) && isempty(function_)
        if length(strfind(tline,'"'))==2
            split=regexp(tline,'\"','split');
            function_ = split(end-1);
            functions(counter)=function_;
        else
            split=regexp(tline,'\"','split');
            functions(counter)=split(end);
            tline=fgetl(fileID);
            split=regexp(tline,'\"','split');
            functions{counter}=[functions{counter} ' ' regexprep(split{1},' ','')];
            function_ = enzymes(counter);
        end
    end
    
    if (~isempty(findstr(tline, '/product="')) || ~isempty(findstr(tline, '/function="'))) && isempty(enzyme)
        if length(strfind(tline,'"'))==2
            split=regexp(tline,'\"','split');
            enzyme = split(end-1);
            enzymes(counter)=enzyme;
        else
            split=regexp(tline,'\"','split');
            enzymes(counter)=split(end);
            tline=fgetl(fileID);
            split=regexp(tline,'\"','split');
            enzymes{counter}=[enzymes{counter} ' ' regexprep(split{1},' ','')];
            enzyme = enzymes(counter);
        end
    end
    
    if (~isempty(findstr(tline, '/EC_number="')) || ~isempty(findstr(tline, '/EC="'))) && isempty(ECnumber)
        split=regexp(tline,'\"','split');
        ECnumber = split(end-1);
        EC(counter)= ECnumber;
    end
    
    if ~isempty(findstr(tline, '/pseudo')) && isempty(pseudo)
        pseudoGenes{counter} = '1';
    end
    
    tline=fgetl(fileID);
end
fclose(fileID);

genes=genes(1:counter);
oldGenes=oldGenes(1:counter);
names=names(1:counter);
EC=EC(1:counter);
enzymes=enzymes(1:counter);
sequences=sequences(1:counter);
pseudoGenes = pseudoGenes(1:counter);
protein_ids = protein_ids(1:counter);
functions = functions(1:counter);
products = products(1:counter);
pseudoGenes( find(cellfun(@isempty, pseudoGenes))) = repmat({'0'},length( find(cellfun(@isempty, pseudoGenes))),1);
pseudoGenes_num = zeros(length(pseudoGenes),1);
pseudoGenes_num( find(strcmp('1',pseudoGenes))) = 1;

info.locus_tag = genes;
info.old_locus_tag = oldGenes;
info.name = names;
info.EC = EC;
info.enzymes = enzymes;
info.sequences = sequences;
info.psedogenes = pseudoGenes_num;
info.functions = functions;
info.products = products;
info.protein_ids = protein_ids;


for i = 1:length(info.old_locus_tag); if isempty(info.old_locus_tag{i}); info.old_locus_tag{i} =''; end; end;
for i = 1:length(info.name); if isempty(info.name{i}); info.name{i} =''; end; end;
for i = 1:length(info.EC); if isempty(info.EC{i}); info.EC{i} =''; end; end;
for i = 1:length(info.enzymes); if isempty(info.enzymes{i}); info.enzymes{i} =''; end; end;
for i = 1:length(info.sequences); if isempty(info.sequences{i}); info.sequences{i} =''; end; end;
for i = 1:length(info.products); if isempty(info.products{i}); info.products{i} =''; end; end;
for i = 1:length(info.functions); if isempty(info.functions{i}); info.functions{i} =''; end; end;
for i = 1:length(info.protein_ids); if isempty(info.protein_ids{i}); info.protein_ids{i} =''; end; end;


end