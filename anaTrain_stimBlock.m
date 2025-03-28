% anaTrain_stimBlock
clear all
close all

%% Generate Dataset

dataFold = 'Y:\Brandon\data';
proj = 'Train_V1Cool_stimBlock';

% projectTbl=getProjectFiles(proj,1,'age','recSite','priorMFlag','priorDescr','duringMFlag','manipDescr','manipDetail');
% for e = 1:height(projectTbl)
% 
%     animal = projectTbl.experimentId{e};
%     unit = projectTbl.unitNr{e};
%     expt = projectTbl.experimentNr{e};
%     exptName = [animal '_u' unit '_' expt];
%     probe = projectTbl.probeId(e);
% 
% %     x(e) = isfile(fullfile(dataFold,'Ephys',animal,exptName,[exptName '_p' num2str(probe) '_MUspkMerge.mat']));
%     if strcmp(exptName,'febk7_u000_009') || strcmp(exptName,'febl7_u000_031')
%         continue
%     end
%     disp(['generating sumStats for ' exptName])
%     [sumStats{e,1}] = anaOri(animal,unit,expt,probe,'MU',dataFold,0,0);
% end
% projectTbl.sumStats = sumStats;
% save(fullfile(dataFold,'dataSets','training',proj(1:end-10),'MU',[proj '_dataSet.mat']),'proj','projectTbl')

load(fullfile(dataFold,'dataSets','training',proj(1:end-10),'MU',[proj '_dataSet.mat']),'projectTbl')
stimBlockTbl = projectTbl;
clear projectTbl

% load(fullfile(dataFold,'dataSets','training',proj(1:end-10),'MU',[proj(1:end-10) '_MU_dataSet.mat']),'projectTbl')
% evalTbl = projectTbl;
% clear projectTbl

%% Analysis

area = 'PSS';
areaIdx = strcmp(stimBlockTbl.recSite,area);
stimBlockTbl = stimBlockTbl(areaIdx,:);

animals = unique(stimBlockTbl.experimentId);
for a = 1:length(animals)
    curAni = animals{a};
    curAniIdx = find(strcmp(stimBlockTbl.experimentId,curAni));
    for b = 1:length(curAniIdx)
        idx = curAniIdx(b);
        blockName{b,a} = [stimBlockTbl.experimentId{idx} '_u' stimBlockTbl.unitNr{idx} '_' stimBlockTbl.experimentNr{idx}];
        blockData{b,a} = stimBlockTbl.sumStats{idx};
    end
end
nBlk = size(blockData,1);
nAni = size(blockData,2);

%% Plot

clrs = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],...
        [0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],...
        [0.6350 0.0780 0.1840]};
figure;hold on
for a = 1:nAni
    for b = 1:nBlk
        if isempty(blockData{b,a})
            nGU_stimBlock(b,a) = nan;
            continue
        end
        nGU_stimBlock(b,a) = sum(blockData{b,a}.goodUnit);
    end
    plot(nGU_stimBlock(:,a),'-o')
end




