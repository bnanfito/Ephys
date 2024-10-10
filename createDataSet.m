%create dataset

clear all
close all

zPath = 'Z:\EphysNew\data';
anaPath = 'Z:\EphysNew\analyzer';
dataPath = 'D:\data';


ageLims = [0 1000];

[coolData] = queryCoolV1(ageLims);
[cntrlData] = queryCntrl(ageLims);

%only take cntrl data that has accompanying cooling data
coolAnimals = unique(coolData.experimentId);
cntrlIdx = ismember(cntrlData.experimentId,coolAnimals);
cntrlData = cntrlData(cntrlIdx,:);
clear cntrlIdx

%% remove Kalatsky experiments
coolKalIdx = strcmp(coolData.module,'BK');
coolKalatskyData = coolData(coolKalIdx,:);
coolData = coolData(~coolKalIdx,:);

cntrlKalIdx = strcmp(cntrlData.module,'BK');
cntrlKalatskyData = cntrlData(cntrlKalIdx,:);
cntrlData = cntrlData(~cntrlKalIdx,:);

%% index experiments

coolLooperVals = [coolData.looperNameCond1 coolData.looperNameCond2 ...
                  coolData.looperNameCond3 coolData.looperNameCond4 ...
                  coolData.looperNameCond5];
cntrlLooperVals = [cntrlData.looperNameCond1 cntrlData.looperNameCond2 ...
                   cntrlData.looperNameCond3 cntrlData.looperNameCond4 ...
                   cntrlData.looperNameCond5];

% contralateral hemifield experiments
stimParams = {'x_size','y_size','x_pos','y_pos'};

hemiCoolIdx = sum(contains(coolLooperVals,'size'),2)>0;
for i = 1:height(coolData)
    fileBase = [coolData.experimentId{i} '_u' coolData.unitNr{i} '_' coolData.experimentNr{i}];
    load(fullfile(anaPath,fileBase(1:5),[fileBase '.analyzer']),'-mat')
    for p = 1:length(stimParams)
        spValCool(i,p) = getparam(stimParams{p},Analyzer);
    end
end
hemiCoolIdx = hemiCoolIdx | spValCool( : , contains(stimParams,'x_size') ) == 50 | ...
                    spValCool( : , contains(stimParams,'y_size') ) == 50;


hemiCntrlIdx = sum(contains(cntrlLooperVals,'size'),2)>0;
for i = 1:height(cntrlData)
    fileBase = [cntrlData.experimentId{i} '_u' cntrlData.unitNr{i} '_' cntrlData.experimentNr{i}];
    load(fullfile(anaPath,fileBase(1:5),[fileBase '.analyzer']),'-mat')
    for p = 1:length(stimParams)
        spValCntrl(i,p) = getparam(stimParams{p},Analyzer);
    end
end
hemiCntrlIdx = hemiCntrlIdx | spValCntrl( : , contains(stimParams,'x_size') ) == 50 | ...
                    spValCntrl( : , contains(stimParams,'y_size') ) == 50;


% contrast experiments
coolContIdx = sum(contains(coolLooperVals,'contrast'),2)>0;
cntrlContIdx = sum(contains(cntrlLooperVals,'contrast'),2)>0;

% ori expts
coolOriIdx = sum(contains(coolLooperVals,{'ori16','ori12'}),2)>0;
cntrlOriIdx = sum(contains(cntrlLooperVals,{'ori16','ori12'}),2)>0;






