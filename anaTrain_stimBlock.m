% anaTrain_stimBlock
clear all
close all

%% Load Project Table

dataFold = 'Y:\Brandon\data';
proj = 'Train_V1Cool_stimBlock';
anaMode = 'MU';

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
%     [sumStats{e,1}] = anaOri(animal,unit,expt,probe,anaMode,dataFold,0,0);
% end
% projectTbl.sumStats = sumStats;
% save(fullfile(dataFold,'dataSets','training',proj(1:end-10),anaMode,[proj '_' anaMode 'dataSet.mat']),'proj','projectTbl')

load(fullfile(dataFold,'dataSets','training',proj(1:end-10),anaMode,[proj '_' anaMode 'dataSet.mat']),'projectTbl')

%% Organize Data

area = 'PSS';
areaIdx = strcmp(projectTbl.recSite,area);
projectTbl = projectTbl(areaIdx,:);

animals = unique(projectTbl.experimentId);
for a = 1:length(animals)
    curAni = animals{a};
    curAniIdx = find(strcmp(projectTbl.experimentId,curAni));
    for b = 1:length(curAniIdx)
        curBlockIdx = curAniIdx(b);
        blockName{b,a} = [projectTbl.experimentId{curBlockIdx} '_u' projectTbl.unitNr{curBlockIdx} '_' projectTbl.experimentNr{curBlockIdx}];
        if ~isempty(projectTbl.sumStats{curBlockIdx})
            blockData{b,a} = projectTbl.sumStats{curBlockIdx}(projectTbl.sumStats{curBlockIdx}.goodUnit,:);
        end
    end
end
nBlk = size(blockData,1);
nAni = size(blockData,2);

%% Plot

clrs = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],...
        [0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],...
        [0.6350 0.0780 0.1840]};
metric = 'rPref';
figure;hold on
for b = 1:nBlk
    gX = 0.6;
    x = b-(gX/2):gX/(nAni+1):b+(gX/2);
    X(b,:) = x(2:end-1);
    for a = 1:nAni
        if isempty(blockData{b,a})
            Y{b,a} = [];
            meanY(b,a) = nan;
        else
            Y{b,a} = blockData{b,a}{:,metric};
            meanY(b,a) = mean(Y{b,a},'omitnan');
            plot(X(b,a),Y{b,a},'.','Color',clrs{a})
        end
    end

end
for a = 1:nAni
    plot(X(:,a),meanY(:,a),'-o','Color',clrs{a})
end
xlabel('stimBlock number')
ylabel(metric)


