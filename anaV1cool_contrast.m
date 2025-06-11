% anaV1cool_contrast
clear
close all

anaMode = 'MU';
proj = 'V1cool_contrast';
area = 'PSS';
dataFold = 'Y:\Brandon\data';
% dataFold = '/Volumes/NielsenHome2/Brandon/data';

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

load(fullfile(dataFold,'dataSets','cooling',proj,anaMode,'V1cool_contrast_MUdataSet.mat'))

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
    coolDat = vertcat(dat.cool{a,:});
    if isempty(cntrlDat) || isempty(coolDat)
        distCntrl{a} = [];
        meanCntrl(a) = nan;
        semCntrl(a) = nan;
        nUcntrl(a) = nan;
        distCool{a} = [];
        meanCool(a) = nan;
        semCool(a) = nan;
        nUcool(a) = nan;
        continue
    end

    goodIdxCntrl = screenUnits(cntrlDat,anaMode);
    distCntrl{a} = cntrlDat(goodIdxCntrl,:).cF;
    meanCntrl(a) = mean(distCntrl{a},'omitnan');
    semCntrl(a) = std(distCntrl{a},'omitnan')/sqrt(length(distCntrl{a}));
    nUcntrl(a) = length(distCntrl{a});

    goodIdxCool = screenUnits(coolDat,anaMode);
    distCool{a} = coolDat(goodIdxCool,:).cF;
    meanCool(a) = mean(distCool{a},'omitnan');
    semCool(a) = std(distCool{a},'omitnan')/sqrt(length(distCool{a}));
    nUcool(a) = length(distCool{a});

end

%% Plot

nUmin = 20;

figure; hold on
plot(ages(nUcntrl>=nUmin),meanCntrl(nUcntrl>=nUmin),'ko','MarkerFaceColor','k')
plot(ages(nUcntrl<nUmin),meanCntrl(nUcntrl<nUmin),'ko','LineWidth',2)
plot(repmat(ages,2,1),meanCntrl+([1;-1]*semCntrl),'k','LineWidth',2)
plot(ages(nUcool>=nUmin),meanCool(nUcool>=nUmin),'co','MarkerFaceColor','c')
plot(ages(nUcool<nUmin),meanCool(nUcool<nUmin),'co','LineWidth',2)
plot(repmat(ages,2,1),meanCool+([1;-1]*semCool),'c','LineWidth',2)
xlabel('age (postnatal day)')
ylabel('C50')
