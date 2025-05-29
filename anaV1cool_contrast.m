% anaV1cool_contrast
clear
close all

anaMode = 'MU';
proj = 'V1cool_contrast';
area = 'PSS';
dataFold = fullfile('Y:\Brandon\data');

%% generate project table

projTbl = getProjectFiles(proj,1,'age','recSite','penNr','priorMFlag','priorDescr',...
                                       'duringMFlag','manipDescr','manipDetail',...
                                       'looperNameCond1','looperNameCond2',...
                                       'looperNameCond3','looperNameCond4',...
                                       'looperNameCond5');
projTbl = projTbl(strcmp(projTbl.recSite,area),:);
animals = unique(projTbl.experimentId);
for a = 1:length(animals)
    aniIdx = strcmp(projTbl.experimentId,animals{a});
    ages(a) = unique(projTbl.age(aniIdx));
end
for e = 1:height(projTbl)
    animal = projTbl.experimentId{e};
    unit = projTbl.unitNr{e};
    expt = projTbl.experimentNr{e};
    probe = projTbl.probeId(e);
    exptName = projTbl.fileBase{e};
    disp(['generating sumStats for ' exptName])
    checkExp(e) = isfile(fullfile(dataFold,'Ephys',animal,exptName,[exptName '_p' num2str(probe) '_MUspkMerge.mat']));
%     [sumStats{e,1}] = anaCon(animal,unit,expt,probe,anaMode,dataFold,0,0);
%     nGU(e,1) = sum(screenUnits(sumStats{e,1},anaMode));
end
needsMerge = {projTbl.fileBase{checkExp==0}};
needsMerge = needsMerge';
