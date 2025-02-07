
clear all
close all

%Settings
animalId = 'febj7';
mergeId = '000015000017000018';
probeId = 1;
mergeName = [animalId '_uMMM_' mergeId];
dataFold = 'F:\Brandon\data';
anaMode = 'SU';

%Load merge info
load(fullfile(dataFold,'Ephys',animalId,mergeName,[mergeName '_mergeInfo.mat']))
nFiles = length(mergeInfo.files);
for f = 1:nFiles
    exptName{f,1} = [animalId '_' mergeInfo.files{f}];

    sumStats{f} = anaOri(animalId,exptName{f}(8:10),exptName{f}(12:14),probeId,anaMode,dataFold,0,0,f);
    
end
