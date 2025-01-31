%anaV1cool_ori2
clear
close all

anaMode = 'MU';
proj = ['V1cool_MU_ori'];
area = 'PSS';

% dataFold = fullfile('Y:\Brandon\data');
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
%     [sumStats{e,1}] = anaOri(animal,unit,expt,probe,anaMode,dataFold,0,0);
%     nGoodUnits(e,1) = sum(sumStats{e}.goodUnit);
% end
% projTbl.sumStats = sumStats;
% projTbl.nGU = nGoodUnits;

% load(['Y:\Brandon\data\dataSets\cooling\' proj '\V1cool_' anaMode '_ori_projectTbl.mat'])
load(['/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/cooling/' proj '/V1cool_' anaMode '_ori_projectTbl.mat'])

coolIdx = projTbl.duringMFlag == 1 & strcmp(projTbl.manipDescr,'Cooling');
cntrlIdx = projTbl.duringMFlag == 0 & projTbl.priorMFlag == 0;
for a = 1:length(animals)
    aniIdx = strcmp(projTbl.experimentId,animals{a});

    cntrlPens = unique(projTbl.penNr(cntrlIdx & aniIdx));
    for p = 1:length(cntrlPens)
        penIdx = projTbl.penNr == cntrlPens(p);
        idx = find(cntrlIdx & aniIdx & penIdx);
        if length(idx) > 1 %if there are repeated expts in a penetration, only take the one with highest yield of good units
            idx = idx( projTbl.nGU(idx) == max(projTbl.nGU(idx)) );
        end
        dat{a,cntrlPens(p),1} = projTbl.sumStats{idx};
    end

    coolPens = unique(projTbl.penNr(coolIdx & aniIdx));
    for p = 1:length(coolPens)
        penIdx = projTbl.penNr == coolPens(p);
        idx = find(coolIdx & aniIdx & penIdx);
        if length(idx) > 1
            idx = idx( projTbl.nGU(idx) == max(projTbl.nGU(idx)) );
        end
        dat{a,coolPens(p),2} = projTbl.sumStats{idx};

    end

end


SI = nan(length(animals),1);
semSI = nan(length(animals),1);
nU = nan(length(animals),1);
meanCntrl = nan(length(animals),1);
semCntrl = nan(length(animals),1);
meanCool = nan(length(animals),1);
semCool = nan(length(animals),1);
for a = 1:length(animals)

    if isempty(dat{a,1,1}) || isempty(dat{a,1,2})
        continue
    end
    goodId = dat{a,1,1}.uID(dat{a,1,1}.goodUnit);
    nU(a,1) = length(goodId);
    cntrlDist = dat{a,1,1}.rPref(ismember(dat{a,1,1}.uID,goodId));
    coolDist = dat{a,1,2}.rPref(ismember(dat{a,1,2}.uID,goodId));
    siDist = (coolDist-cntrlDist)./cntrlDist;

    SI(a) = mean(siDist,'omitnan');
    semSI(a) = std(siDist,'omitnan')/sqrt(length(siDist));
    meanCntrl(a) = mean(cntrlDist,'omitnan');
    semCntrl(a) = std(cntrlDist,'omitnan')/sqrt(length(cntrlDist));
    meanCool(a) = mean(coolDist,'omitnan');
    semCool(a) = std(coolDist,'omitnan')/sqrt(length(coolDist));

end

uMin = 20;

figure; hold on
plot(ages(nU>=uMin),SI(nU>=uMin),'k.','MarkerSize',20)
plot(ages(nU<uMin),SI(nU<uMin),'ko','MarkerSize',7,'LineWidth',2)
plot(repmat(ages,2,1),SI'+([1;-1]*semSI'),'k','LineWidth',2)
xlabel('age (postnatal day)')
ylabel('delta rPref: (cool- cntrl)/cntrl')

figure; hold on
plot(ages(nU>=uMin),meanCntrl(nU>=uMin),'k.','MarkerSize',20)
plot(ages(nU<uMin),meanCntrl(nU<uMin),'ko','MarkerSize',7,'LineWidth',2)
plot(repmat(ages,2,1),meanCntrl'+([1;-1]*semCntrl'),'k','LineWidth',2)
plot(ages(nU>=uMin),meanCool(nU>=uMin),'c.','MarkerSize',20)
plot(ages(nU<uMin),meanCool(nU<uMin),'co','MarkerSize',7,'LineWidth',2)
plot(repmat(ages,2,1),meanCool'+([1;-1]*semCool'),'c','LineWidth',2)
xlabel('age (postnatal day)')
ylabel('rPref (Hz)')

