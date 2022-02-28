function downloadVHMModels

% [n,s] = xlsread('microbes');
% modelsID = s(2:end,29);
load('modelsID')
cd('C:\Users\notebook\Desktop\Trabajo\Gut\AllModels\MAT')
for i = 1:length(modelsID)
    if exist(['C:\Users\notebook\Desktop\Trabajo\Gut\AllModels\MAT', filesep, modelsID{i}, '.mat'], 'file') ~= 2
        Direction = ['https://webdav-r3lab.uni.lu/public/msp/AGORA/mat/' modelsID{i} '.mat'];
        [status_curl, result_curl] = system(['curl --max-time 15 -s -k -L --head ', Direction]);
        
        % check if the URL exists
        if status_curl == 0 && ~isempty(strfind(result_curl, '200 OK'))
            status_curlDownload = system(['curl --max-time 6000 -O -L ', Direction]);
            
            if status_curlDownload == 0
                fprintf(' + Downloaded:      %s\n', Direction);
            end
        else
            fprintf(' > The URL %s cannot be reached.\n', Direction);
        end
    end
end

cd('C:\Users\notebook\Desktop\Trabajo\Gut')

end