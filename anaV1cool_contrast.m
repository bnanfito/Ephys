% anaV1cool_contrast
clear
close all

anaMode = 'MU';
proj = 'V1cool_contrast';
area = 'PSS';
% dataFold = 'Y:\Brandon\data';
dataFold = '/Volumes/NielsenHome2/Brandon/data';

%% generate project table

% projTbl = getProjectFiles(proj,1,'age','recSite','penNr','priorMFlag','priorDescr',...
%                                        'duringMFlag','manipDescr','manipDetail',...
%                                        'looperNameCond1','looperNameCond2',...
%                                        'looperNameCond3','looperNameCond4',...
%                                        'looperNameCond5');
% projTbl = projTbl(strcmp(projTbl.recSite,area),:);
% animals = unique(projTbl.experimentId);
% for a = 1:length(animals)
%     aniIdx = strcmp(projTbl.experimentId,animals{a});
%     ages(a) = unique(projTbl.age(aniIdx));
% end
% for e = 1:height(projTbl)
%     animal = projTbl.experimentId{e};
%     unit = projTbl.unitNr{e};
%     expt = projTbl.experimentNr{e};
%     probe = projTbl.probeId(e);
%     exptName = projTbl.fileBase{e};
%     disp(['generating sumStats for ' exptName])
%     [sumStats{e,1}] = anaCon(animal,unit,expt,probe,anaMode,dataFold,0,0);
%     nGU(e,1) = sum(screenUnits(sumStats{e,1},anaMode));
% end
% projTbl.sumStats = sumStats;
% projTbl.nGU = nGU;

load(fullfile(dataFold,'dataSets','cooling','V1cool_contrast','V1cool_contrast_MUdataSet.mat'))

%% Organize Data

coolIdx = projTbl.duringMFlag == 1 & strcmp(projTbl.manipDescr,'Cooling');
cntrlIdx = projTbl.duringMFlag == 0 & projTbl.priorMFlag == 0;
for a = 1:length(animals)
    aniIdx = strcmp(projTbl.experimentId,animals{a});

    %if there are repeated expts in a penetration, only take the one with highest yield of good units
    cntrlPens = unique(projTbl.penNr(cntrlIdx & aniIdx));
    for p = 1:length(cntrlPens)
        penIdx = projTbl.penNr == cntrlPens(p);
        idx = find(cntrlIdx & aniIdx & penIdx);
        if length(idx) > 1 
            idx = idx( projTbl.nGU(idx) == max(projTbl.nGU(idx)) );
        end
        dat.cntrl{a,cntrlPens(p)} = projTbl.sumStats{idx};
    end

    coolPens = unique(projTbl.penNr(coolIdx & aniIdx));
    for p = 1:length(coolPens)
        penIdx = projTbl.penNr == coolPens(p);
        idx = find(coolIdx & aniIdx & penIdx);
        if length(idx) > 1
            idx = idx( projTbl.nGU(idx) == max(projTbl.nGU(idx)) );
        end
        dat.cool{a,coolPens(p)} = projTbl.sumStats{idx};

    end

end

%% Calculate Metrics

for a = 1:length(animals)

    cntrlDat = vertcat(dat.cntrl{a,:});
    if isempty(cntrlDat)
        distCntrl{a} = [];
        meanCntrl(a) = nan;
        stdCntrl(a) = nan;
        continue
    end
    goodIdxCntrl = screenUnits(cntrlDat,anaMode);
    distCntrl{a} = cntrlDat(goodIdxCntrl,:).cF;
    meanCntrl(a) = mean(distCntrl{a},'omitnan');
    stdCntrl(a) = std(distCntrl{a},'omitnan');

    coolDat = vertcat(dat.cool{a,:});
    if isempty(coolDat)
        distCool{a} = [];
        meanCool(a) = nan;
        stdCool(a) = nan;
        continue
    end
    goodIdxCool = screenUnits(coolDat,anaMode);
    distCool{a} = coolDat(goodIdxCool,:).cF;
    meanCool(a) = mean(distCool{a},'omitnan');
    stdCool(a) = std(distCool{a},'omitnan');

end

%% Plot

figure; hold on
for a = 1:length(animals)

    plot(ages(a),meanCntrl(a),'ko')
    plot(ages(a),meanCool(a),'co')

end
