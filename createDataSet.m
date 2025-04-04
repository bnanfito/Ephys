%create dataset

% clear all
% close all

zPath = 'Z:\EphysNew\data';
anaPath = 'Z:\EphysNew\analyzer';
procSpksPath = 'Z:\EphysNew\processedSpikes';
% dataPath = 'F:\Brandon\data';
dest = 'D:\data\dataSets\DSdev\Ephys';

for f = 1:height(projList)
    animal = projList.experimentId{f};
    unit = projList.unitNr{f};
    expt = projList.experimentNr{f};
    probe = projList.probeId(f);
    exptName = [animal '_u' unit '_' expt];

    if ~isfolder(fullfile(dest,animal,exptName))
        mkdir(fullfile(dest,animal,exptName))
    end

    copyfile(fullfile(anaPath,animal,[exptName '.analyzer']),fullfile(dest,animal,exptName,[exptName '.analyzer']))
    copyfile(fullfile(procSpksPath,animal,exptName,[exptName '_id.mat']),fullfile(dest,animal,exptName,[exptName '_id.mat']))
    copyfile(fullfile(procSpksPath,animal,exptName,[exptName '_trialInfo.mat']),fullfile(dest,animal,exptName,[exptName '_trialInfo.mat']))
    copyfile(fullfile(procSpksPath,animal,exptName,[exptName '_p' num2str(probe) '_spkSort.mat']),fullfile(dest,animal,exptName,[exptName '_p' num2str(probe) '_spkSort.mat']))

end









% ageLims = [0 1000];
% 
% [coolData] = queryCoolV1(ageLims);
% [cntrlData] = queryCntrl(ageLims);
% 
% % %only take cntrl data that has accompanying cooling data
% % coolAnimals = unique(coolData.experimentId);
% % cntrlIdx = ismember(cntrlData.experimentId,coolAnimals);
% % cntrlData = cntrlData(cntrlIdx,:);
% % clear cntrlIdx
% % 
% % %% remove Kalatsky experiments
% % coolKalIdx = strcmp(coolData.module,'BK');
% % coolKalatskyData = coolData(coolKalIdx,:);
% % coolData = coolData(~coolKalIdx,:);
% % 
% % cntrlKalIdx = strcmp(cntrlData.module,'BK');
% % cntrlKalatskyData = cntrlData(cntrlKalIdx,:);
% % cntrlData = cntrlData(~cntrlKalIdx,:);
% % 
% % %% index experiments
% % 
% % coolLooperVals = [coolData.looperNameCond1 coolData.looperNameCond2 ...
% %                   coolData.looperNameCond3 coolData.looperNameCond4 ...
% %                   coolData.looperNameCond5];
% % cntrlLooperVals = [cntrlData.looperNameCond1 cntrlData.looperNameCond2 ...
% %                    cntrlData.looperNameCond3 cntrlData.looperNameCond4 ...
% %                    cntrlData.looperNameCond5];
% % 
% % % contralateral hemifield experiments
% % stimParams = {'x_size','y_size','x_pos','y_pos'};
% % 
% % hemiCoolIdx = sum(contains(coolLooperVals,'size'),2)>0;
% % for i = 1:height(coolData)
% %     fileBase = [coolData.experimentId{i} '_u' coolData.unitNr{i} '_' coolData.experimentNr{i}];
% %     load(fullfile(anaPath,fileBase(1:5),[fileBase '.analyzer']),'-mat')
% %     for p = 1:length(stimParams)
% %         spValCool(i,p) = getparam(stimParams{p},Analyzer);
% %     end
% % end
% % hemiCoolIdx = hemiCoolIdx | spValCool( : , contains(stimParams,'x_size') ) == 50 | ...
% %                     spValCool( : , contains(stimParams,'y_size') ) == 50;
% % 
% % 
% % hemiCntrlIdx = sum(contains(cntrlLooperVals,'size'),2)>0;
% % for i = 1:height(cntrlData)
% %     fileBase = [cntrlData.experimentId{i} '_u' cntrlData.unitNr{i} '_' cntrlData.experimentNr{i}];
% %     load(fullfile(anaPath,fileBase(1:5),[fileBase '.analyzer']),'-mat')
% %     for p = 1:length(stimParams)
% %         spValCntrl(i,p) = getparam(stimParams{p},Analyzer);
% %     end
% % end
% % hemiCntrlIdx = hemiCntrlIdx | spValCntrl( : , contains(stimParams,'x_size') ) == 50 | ...
% %                     spValCntrl( : , contains(stimParams,'y_size') ) == 50;
% % 
% % 
% % % contrast experiments
% % coolContIdx = sum(contains(coolLooperVals,'contrast'),2)>0;
% % cntrlContIdx = sum(contains(cntrlLooperVals,'contrast'),2)>0;
% % 
% % % ori expts
% % coolOriIdx = sum(contains(coolLooperVals,{'ori16','ori12'}),2)>0;
% % cntrlOriIdx = sum(contains(cntrlLooperVals,{'ori16','ori12'}),2)>0;
% 
% area = 'PSS';
% looperVars = [cntrlData.looperNameCond1 cntrlData.looperNameCond2 cntrlData.looperNameCond3 cntrlData.looperNameCond4 cntrlData.looperNameCond5];
% areaIdx = strcmp(cntrlData.recSite,area);
% contIdx = sum(contains(looperVars,'contrast'),2)>0;
% % take only contrast expt & recordings from 'area'
% cntrlData = cntrlData(areaIdx & contIdx,:);
% 
% for i = 1:height(cntrlData)
%     
%     exptId = cntrlData.experimentId{i};
%     exptName = [exptId '_u' cntrlData.unitNr{i} '_' cntrlData.experimentNr{i}];
%     load(fullfile(anaPath,exptId,[exptName '.analyzer']),'-mat')
%     xSize(i) = getparam('x_size',Analyzer);
%     ySize(i) = getparam('y_size',Analyzer);
%     noMaskBit(i) = strcmp(getparam('mask_type',Analyzer),'none');
%     clear Analyzer
% end
% ffIdx = xSize>=100 & ySize>=100 & noMaskBit;
% % take only full field experiments
% cntrlData = cntrlData(ffIdx,:);
% 
% animals = unique(cntrlData.experimentId);
% for a = 1:length(animals)
%     ages(a) = unique(cntrlData.age(strcmp(cntrlData.experimentId,animals{a})));
% end
% 
% figure; hold on
% histogram(ages,[0:1:1000])
% 
% fileName = fullfile(dataPath,'dataSets','cntrl','cntrlContrast.xlsx');
% writetable(cntrlData,fileName,'Sheet',1,'Range','A1')
% 
% fileName = fullfile(dataPath,'dataSets','cntrl','cntrlContrast.mat');
% save(fileName,'cntrlData','animals','ages')




