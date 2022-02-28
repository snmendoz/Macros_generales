function downloadBIGGDatabase
if exist('D:\Dropbox\Databases\BIGG\logBIGGdownload','file')==2
%     delete('logBIGGdownload')
else
    d = datetime;
    nameDiary = ['logBIGGdownload_' regexprep(datestr(d),{':',' ','-'},{'_','_','_'})];
    diary(nameDiary);
    fprintf('Executed on:      %s\n', d)
    str = urlread('http://bigg.ucsd.edu/api/v2/models');
    pos = strfind(str,'"bigg_id": ');
    ids =cell(length(pos),1);
    for i = 1:length(pos)
        if i<length(pos)
            substr = str(pos(i)+1:pos(i+1)-1);
        else
            substr = str(pos(i)+1:end);
        end
        pos2 = strfind(substr,'"');
        ids{i} = substr(pos2(2)+1:pos2(3)-1);
    end
    fi = fopen('modelNames.txt', 'w+');
    for i = 1:length(ids);
        fprintf(fi, [ ids{i} '\n']);
    end
    fclose(fi);
    
    printLevel =1;
    if printLevel == 0
        curlSilence = '-s';
    else
        curlSilence = '';
    end
    
    modelArr = strcat('http://bigg.ucsd.edu/static/models/',ids,'.mat');
    
    cd('D:\Dropbox\Databases\BIGG');
    if exist('D:\Dropbox\Databases\BIGG\succesfullyDownloaded.mat','file')==2
        load('succesfullyDownloaded.mat'); 
        if length(modelArr)>length(succesfullyDownloaded)
        succesfullyDownloaded = [succesfullyDownloaded zeros(1,length(modelArr)-length(succesfullyDownloaded))];
        end
    else
        succesfullyDownloaded = zeros(size(pos));
    end
    for i = 1:length(modelArr)
        if succesfullyDownloaded(i) && exist([ids{i} '.mat'],'file')==2
        else
            
            [status_curl, result_curl] = system(['curl --max-time 15 -s -k -L --head ', modelArr{i}]);
            
            % check if the URL exists
            if status_curl == 0 && ~isempty(strfind(result_curl, '200 OK'))
                status_curlDownload = -1;
                n_attempts = 0;
                
                while status_curlDownload ~=0 && n_attempts < 5
                    n_attempts = n_attempts + 1;
                    status_curlDownload = system(['curl ', curlSilence, ' --max-time 300 -O -L ', modelArr{i}]);
                    
                    if printLevel > 0 && status_curlDownload == 0
                        fprintf(' + Downloaded:      %s\n', modelArr{i});
                        succesfullyDownloaded(i) = 1;
                    end
                end
                
            else
                fprintf(' > The URL %s cannot be reached.\n', modelArr{i});
                
            end
        end
    end
    save('succesfullyDownloaded', 'succesfullyDownloaded');
    diary off
end