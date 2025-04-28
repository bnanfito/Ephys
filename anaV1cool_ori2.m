%anaV1cool_ori2
clear
close all

anaMode = 'MU';
proj = 'V1cool_ori';
area = 'PSS';

%% Generate Project Table

dataFold = fullfile('Y:\Brandon\data');
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
%     nGU(e,1) = sum(screenUnits(sumStats{e,1},anaMode));
% end
% projTbl.sumStats = sumStats;
% projTbl.nGU = nGU;

load(fullfile(dataFold,'dataSets','cooling',proj,anaMode,[proj '_' anaMode 'dataSet.mat']))
% load(['/Volumes/NielsenHome2/Brandon/data/dataSets/cooling/' proj '/V1cool_' anaMode '_ori_projectTbl.mat'])
% load(['/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/cooling/' proj '/V1cool_' anaMode '_ori_projectTbl.mat'])

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

SI = nan(length(animals),1);
semSI = nan(length(animals),1);
nU = nan(length(animals),1);
meanCntrl = nan(length(animals),1);
semCntrl = nan(length(animals),1);
meanCool = nan(length(animals),1);
semCool = nan(length(animals),1);
for a = 1:length(animals)

    if strcmp(anaMode,'MU') && ~(isempty(dat.cntrl{a,1}) || isempty(dat.cool{a,1}))

%         % Suppresion index
%         goodId = dat.cntrl{a,1}.uID(screenUnits(dat.cntrl{a,1},anaMode)); %unit ids that pass inclusion criteria in cntrl condition
%         nU(a,1) = length(goodId);
%         cntrlDist = dat.cntrl{a,1}.rPref(ismember(dat.cntrl{a,1}.uID,goodId));
%         coolDist = dat.cool{a,1}.rPref(ismember(dat.cool{a,1}.uID,goodId));
%         siDist = (coolDist)./(cntrlDist);
%         SI(a) = mean(siDist,'omitnan');
%         semSI(a) = std(siDist,'omitnan')/sqrt(length(siDist));
%         meanCntrl(a) = mean(cntrlDist,'omitnan');
%         semCntrl(a) = std(cntrlDist,'omitnan')/sqrt(length(cntrlDist));
%         meanCool(a) = mean(coolDist,'omitnan');
%         semCool(a) = std(coolDist,'omitnan')/sqrt(length(coolDist));

        % Delta Ldir
        goodId_cntrl = dat.cntrl{a,1}.uID(screenUnits(dat.cntrl{a,1},anaMode));
        goodId_cool = dat.cool{a,1}.uID(screenUnits(dat.cool{a,1},anaMode));
        goodId = goodId_cntrl(ismember(goodId_cntrl,goodId_cool));
        nU(a,1) = length(goodId);
        nU_cntrl(a,1) = length(goodId_cntrl);
        nU_cool(a,1) = length(goodId_cool);
        cntrlDist = dat.cntrl{a,1}.ldr(screenUnits(dat.cntrl{a,1},anaMode));
        coolDist = dat.cool{a,1}.ldr(screenUnits(dat.cool{a,1},anaMode));
        dLdrDist = dat.cool{a,1}.ldr(ismember(dat.cool{a,1}.uID,goodId))-dat.cntrl{a,1}.ldr(ismember(dat.cntrl{a,1}.uID,goodId));
        dLDR(a) = mean(dLdrDist,'omitnan');
        sem_dLDR(a) = std(dLdrDist,'omitnan')/sqrt(length(dLdrDist));
        meanCntrl(a) = mean(cntrlDist,'omitnan');
        semCntrl(a) = std(cntrlDist,'omitnan')/sqrt(length(cntrlDist));
        meanCool(a) = mean(coolDist,'omitnan');
        semCool(a) = std(coolDist,'omitnan')/sqrt(length(coolDist));

    end



end

%% Plot

uMin = 20;

if strcmp(anaMode,'MU')

% figure; hold on
% plot(ages(nU>=uMin),SI(nU>=uMin),'k.','MarkerSize',20)
% plot(ages(nU<uMin),SI(nU<uMin),'ko','MarkerSize',7,'LineWidth',2)
% plot(repmat(ages,2,1),SI'+([1;-1]*semSI'),'k','LineWidth',2)
% xlabel('age (postnatal day)')
% ylabel('SI: Rcool/Rcntrl')
% 
% figure; hold on
% plot(ages(nU>=uMin),meanCntrl(nU>=uMin),'k.','MarkerSize',20)
% plot(ages(nU<uMin),meanCntrl(nU<uMin),'ko','MarkerSize',7,'LineWidth',2)
% plot(repmat(ages,2,1),meanCntrl'+([1;-1]*semCntrl'),'k','LineWidth',2)
% plot(ages(nU>=uMin),meanCool(nU>=uMin),'c.','MarkerSize',20)
% plot(ages(nU<uMin),meanCool(nU<uMin),'co','MarkerSize',7,'LineWidth',2)
% plot(repmat(ages,2,1),meanCool'+([1;-1]*semCool'),'c','LineWidth',2)
% xlabel('age (postnatal day)')
% ylabel('rPref (Hz)')

figure; hold on
plot(ages(nU>=uMin),dLDR(nU>=uMin),'k.','MarkerSize',20)
plot(ages(nU<uMin),dLDR(nU<uMin),'ko','MarkerSize',7,'LineWidth',2)
plot(repmat(ages,2,1),dLDR+([1;-1]*sem_dLDR),'k','LineWidth',2)
yline(0,'k--')
xlabel('age (postnatal day)')
ylabel('delta Ldir')

figure; hold on
plot(ages(nU_cntrl>=uMin),meanCntrl(nU_cntrl>=uMin),'k.','MarkerSize',20)
plot(ages(nU_cntrl<uMin),meanCntrl(nU_cntrl<uMin),'ko','MarkerSize',7,'LineWidth',2)
plot(repmat(ages,2,1),meanCntrl'+([1;-1]*semCntrl'),'k','LineWidth',2)
plot(ages(nU_cool>=uMin),meanCool(nU_cool>=uMin),'c.','MarkerSize',20)
plot(ages(nU_cool<uMin),meanCool(nU_cool<uMin),'co','MarkerSize',7,'LineWidth',2)
plot(repmat(ages,2,1),meanCool'+([1;-1]*semCool'),'c','LineWidth',2)
xlabel('age (postnatal day)')
ylabel('Ldir')

end





