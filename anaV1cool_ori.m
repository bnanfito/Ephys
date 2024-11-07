%anaV1Cool_ori
clear
close all

% anaMode = 'MU';
% proj = ['V1cool_' anaMode '_ori'];
% area = 'PSS';
% dataFold = fullfile('Y:\Brandon\data');
% projTbl = getProjectFiles(proj,1,'age','recSite','penNr','priorMFlag','priorDescr',...
%                                        'duringMFlag','manipDescr','manipDetail',...
%                                        'looperNameCond1','looperNameCond2',...
%                                        'looperNameCond3','looperNameCond4',...
%                                        'looperNameCond5');
% projTbl = projTbl(strcmp(projTbl.recSite,area),:);
% animals = unique(projTbl.experimentId);
% for a = 1:length(animals)
%     aniIdx(:,a) = strcmp(projTbl.experimentId,animals{a});
%     ages(a) = unique(projTbl.age(aniIdx(:,a)));
% end
% for e = 1:height(projTbl)
%     animal = projTbl.experimentId{e};
%     unit = projTbl.unitNr{e};
%     expt = projTbl.experimentNr{e};
%     probe = projTbl.probeId(e);
%     exptName = projTbl.fileBase{e};
%     disp(['generating sumStats for ' exptName])
%     [sumStats{e,1}] = anaOri(animal,unit,expt,probe,anaMode,dataFold,0,0);
%     nGoodUnits(e,1) = sum(sumStats{e}.goodUnit);
% end
% projTbl.sumStats = sumStats;
% projTbl.nGU = nGoodUnits;

load('Y:\Brandon\data\dataSets\cooling\V1cool_MU_ori\V1cool_MU_ori_projectTbl.mat')

coolIdx = projTbl.duringMFlag==1 & strcmp(projTbl.manipDescr,'Cooling') & ...
          strcmp(projTbl.manipDetail,'V1');
newTbl = [];
for a = 1:length(animals)
    %cntrl

    %cool

end

%% Plot

figure;
histogram(ages,[20:2:100])

