
%test - anaOri

function [sumStats] = anaOriTest(animal,unit,expt,probe,anaMode,dataFold,plt,svePlt)
%% settings
plr = 0;
alignBit = 1;
stimMode = 'mono c hemi';
visTest = 'ranksum';
alpha = 0.01;

%% load 

physDir = fullfile(dataFold,'Ephys');
figDir = fullfile(dataFold,'Figures');
exptDir = fullfile(physDir,animal,exptName);
load(fullfile(exptDir, [exptName '_id.mat'] ))
load(fullfile(exptDir, [exptName '_trialInfo.mat'] ))
load(fullfile(exptDir, [exptName '.analyzer']),'-mat')

%% compute

if

elseif

end


end