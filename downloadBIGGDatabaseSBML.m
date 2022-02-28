function downloadBIGGDatabaseSBML
cd('D:\Dropbox\Databases\BIGG');
if exist('D:\Dropbox\Databases\BIGG\logBIGGSBMLdownload','file')==2
    delete('logBIGGSBMLdownload')
else
    diary logBIGGSBMLdownload
    d = datetime;
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
    
    modelArr = strcat('http://bigg.ucsd.edu/static/models/',ids,'.xml');
    
    if exist('D:\Dropbox\Databases\BIGG\succesfullyDownloaded2.mat','file')==2
        load('succesfullyDownloaded2.mat')
    else
        succesfullyDownloaded2 = zeros(size(pos));
    end
    for i = 1:length(modelArr)
        if succesfullyDownloaded2(i) && exist([ids{i} '.xml'],'file')==2
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
                        succesfullyDownloaded2(i) = 1;
                    end
                end
            else
                fprintf(' > The URL %s cannot be reached.\n', modelArr{i});
                
            end
        end
    end
    save('succesfullyDownloaded2', 'succesfullyDownloaded2');
    diary off
end