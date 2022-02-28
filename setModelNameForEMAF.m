function setModelNameForEMAF(modelName, folder)

cd(folder)

copyfile('D:\Dropbox\emaf\tutorial\runMedia3.py', folder);
copyfile('D:\Dropbox\emaf\tutorial\pushRunMedia3.py', folder);

fi2 = fopen('runMediaAux.py','w+');
fi = fopen('runMedia3.py','r+');

tline = fgetl(fi);
while ischar(tline)
    
    if length(tline)>13 && strcmp(tline(1:13),'modelFile = ''')
        fprintf(fi2, ['modelFile = ''' modelName '''\n']);
    else
        fprintf(fi2,[regexprep(tline,{'\\','%'},{'\\\\','%%'}) '\n']);
    end
    tline = fgetl(fi);

end
fclose(fi);
fclose(fi2);
pause(2)
delete('runMedia3.py');
movefile('runMediaAux.py', 'runMedia3.py');

pause(2)
fi3 = fopen('pushMunMediaAux.py','w+');
fi4 = fopen('pushRunMedia3.py','r+');

tline = fgetl(fi4);
while ischar(tline)
    
    if length(tline)>13 && strcmp(tline(1:13),'modelFile = ''')
        fprintf(fi3, ['modelFile = ''' modelName '''\n']);
    else
        fprintf(fi3,[regexprep(tline,{'\\','%'},{'\\\\','%%'}) '\n']);
    end
    tline = fgetl(fi4);

end
fclose(fi4);
fclose(fi3);
pause(2)
delete('pushRunMedia3.py');
movefile('pushMunMediaAux.py', 'pushRunMedia3.py');

end

