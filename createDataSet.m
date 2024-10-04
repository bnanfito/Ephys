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

%% remove Kalatsky experiments
coolKalIdx = strcmp(coolData.module,'BK');
coolKalatskyData = coolData(coolKalIdx,:);
coolData = coolData(~coolKalIdx,:);

cntrlKalIdx = strcmp(cntrlData.module,'BK');
cntrlKalatskyData = cntrlData(cntrlKalIdx,:);
cntrlData = cntrlData(~cntrlKalIdx,:);

%% remove fullfield experiments
% include experiments with size as a looped parameter (should have both ff and hemifield trials)
hemiCoolIdx = contains(coolData.looperNameCond1,'size') | ...
           contains(coolData.looperNameCond2,'size') | ...
           contains(coolData.looperNameCond3,'size') | ...
           contains(coolData.looperNameCond4,'size') | ...
           contains(coolData.looperNameCond5,'size');
stimParams = {'x_size','y_size','x_pos','y_pos'};
for i = 1:height(coolData)
    fileBase = [coolData.experimentId{i} '_u' coolData.unitNr{i} '_' coolData.experimentNr{i}];
    load(fullfile(anaPath,fileBase(1:5),[fileBase '.analyzer']),'-mat')
    spValCool(i,1) = getparam(stimParams{1},Analyzer);
    spValCool(i,2) = getparam(stimParams{2},Analyzer);
    spValCool(i,3) = getparam(stimParams{3},Analyzer);
    spValCool(i,4) = getparam(stimParams{4},Analyzer);
end
hemiCoolIdx = hemiCoolIdx | spValCool( : , contains(stimParams,'x_size') ) == 50 | ...
                    spValCool( : , contains(stimParams,'y_size') ) == 50;
ffCoolData = coolData(~hemiCoolIdx,:);
coolData = coolData(hemiCoolIdx,:);

hemiCntrlIdx = contains(cntrlData.looperNameCond1,'size') | ...
           contains(cntrlData.looperNameCond2,'size') | ...
           contains(cntrlData.looperNameCond3,'size') | ...
           contains(cntrlData.looperNameCond4,'size') | ...
           contains(cntrlData.looperNameCond5,'size');
for i = 1:height(cntrlData)
    fileBase = [cntrlData.experimentId{i} '_u' cntrlData.unitNr{i} '_' cntrlData.experimentNr{i}];
    load(fullfile(anaPath,fileBase(1:5),[fileBase '.analyzer']),'-mat')
    spValCntrl(i,1) = getparam(stimParams{1},Analyzer);
    spValCntrl(i,2) = getparam(stimParams{2},Analyzer);
    spValCntrl(i,3) = getparam(stimParams{3},Analyzer);
    spValCntrl(i,4) = getparam(stimParams{4},Analyzer);
end
hemiCntrlIdx = hemiCntrlIdx | spValCntrl( : , contains(stimParams,'x_size') ) == 50 | ...
                    spValCntrl( : , contains(stimParams,'y_size') ) == 50;
ffCntrlData = cntrlData(~hemiCntrlIdx,:);
cntrlData = cntrlData(hemiCntrlIdx,:);


%% split datasets by domval

coolLooperVals = [coolData.looperNameCond1 coolData.looperNameCond2 ...
                  coolData.looperNameCond3 coolData.looperNameCond4 ...
                  coolData.looperNameCond5];

cntrlLooperVals = [cntrlData.looperNameCond1 cntrlData.looperNameCond2 ...
                   cntrlData.looperNameCond3 cntrlData.looperNameCond4 ...
                   cntrlData.looperNameCond5];

%contrast dataset
cntrlContIdx = sum(contains(cntrlLooperVals,'contrast'),2)>0;
cntrlContData = cntrlData(cntrlContIdx,:);

coolContIdx = sum(contains(coolLooperVals,'contrast'),2)>0;
coolContData = coolData(coolContIdx,:);




