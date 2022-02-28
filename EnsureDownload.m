function EnsureDownload
diary vhmModels2.txt
for i = 1:30
    downloadVHMModels
end
diary off
buildVHMDatabase

end