function newGrRule = replaceGeneInGrRule(grRule,gene, newGene)

condition = ['^' gene '$|^' gene '(?=\ )|(?<=\ )' gene '$|(?<=\ )'  gene '(?=\ )|(' gene '(?=\ )|(?<=\ )' gene ')'];
newGrRule = regexprep(grRule,condition,newGene);

end