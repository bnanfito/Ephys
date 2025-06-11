%anaV1cool_ori2
clear
close all

anaMode = 'MU';
proj = 'V1cool_ori';
area = 'PSS';
dataFold = 'Y:\Brandon\data';
% dataFold = '/Volumes/NielsenHome2/Brandon/data';

%% Generate Project Table

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
% load(fullfile(dataFold,'dataSets','cooling','V1cool_ori','MU','V1cool_ori_MUprojTbl.mat'))
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
        datCntrl = projTbl.sumStats{idx};
        goodId = screenUnits(datCntrl,anaMode);
%         datCntrl = datCntrl(goodId,:);
        dat.cntrl{a,cntrlPens(p)} = datCntrl;
        goodIdCntrl{a,cntrlPens(p)} = goodId;

    end

    coolPens = unique(projTbl.penNr(coolIdx & aniIdx));
    for p = 1:length(coolPens)
        penIdx = projTbl.penNr == coolPens(p);
        idx = find(coolIdx & aniIdx & penIdx);
        if length(idx) > 1
            idx = idx( projTbl.nGU(idx) == max(projTbl.nGU(idx)) );
        end
        datCool = projTbl.sumStats{idx};
        goodId = screenUnits(datCool,anaMode);
%         datCool = datCool(goodId,:);
        dat.cool{a,coolPens(p)} = datCool;
        goodIdCool{a,coolPens(p)} = goodId;

    end

end

%% Calculate Metrics

metric = 'rPref';

SI = nan(length(animals),1);
semSI = nan(length(animals),1);
meanCntrl = nan(length(animals),1);
semCntrl = nan(length(animals),1);
meanCool = nan(length(animals),1);
semCool = nan(length(animals),1);
for a = 1:length(animals)

    if strcmp(anaMode,'MU') && ~(isempty(dat.cntrl{a,1}) || isempty(dat.cool{a,1}))

        % Distributions
        if strcmp(metric,'response var')
            cntrlDist = squeeze(mean(std(cat(3,dat.cntrl{a,1}.response{:}),'omitnan'),'omitnan'));
            coolDist = squeeze(mean(std(cat(3,dat.cool{a,1}.response{:}),'omitnan'),'omitnan'));
        else
            cntrlDist = dat.cntrl{a,1}{:,strcmp(dat.cntrl{a,1}.Properties.VariableNames,metric)}(goodIdCntrl{a,1});
            coolDist = dat.cool{a,1}{:,strcmp(dat.cool{a,1}.Properties.VariableNames,metric)}(goodIdCntrl{a,1});
        end
        meanCntrl(a) = mean(cntrlDist,'omitnan');
        semCntrl(a) = std(cntrlDist,'omitnan')/sqrt(length(cntrlDist));
        meanCool(a) = mean(coolDist,'omitnan');
        semCool(a) = std(coolDist,'omitnan')/sqrt(length(coolDist));
        nUcntrl(a,1) = length(cntrlDist);
        nUcool(a,1) = length(coolDist);

        if strcmp(metric,'rPref')
            % Suppression Index
            siDist{a} = (coolDist-cntrlDist)./(coolDist+cntrlDist);
            siDist_age{a} = repmat(ages(a),length(siDist{a}),1);
            SI(a) = mean(siDist{a},'omitnan');
            semSI(a) = std(siDist{a},'omitnan')/sqrt(length(siDist{a}));
        else
            % Delta Ldir
            goodId_cntrl = dat.cntrl{a,1}.uID(screenUnits(dat.cntrl{a,1},anaMode));
            goodId_cool = dat.cool{a,1}.uID(screenUnits(dat.cool{a,1},anaMode));
            goodId = goodId_cntrl(ismember(goodId_cntrl,goodId_cool));
            nU(a,1) = length(goodId);
            dMetDist = coolDist-cntrlDist;
            dMet(a) = mean(dMetDist,'omitnan');
            sem_dMet(a) = std(dMetDist,'omitnan')/sqrt(length(dMetDist));
        end

    end



end
nID = ~isnan(SI);

%% Plot

uMin = 20;

if strcmp(anaMode,'MU')

% Plot distributions
figure; hold on
plot(ages(nUcntrl>=uMin),meanCntrl(nUcntrl>=uMin),'k.','MarkerSize',20)
plot(ages(nUcntrl<uMin),meanCntrl(nUcntrl<uMin),'ko','MarkerSize',7,'LineWidth',2)
plot(repmat(ages,2,1),meanCntrl'+([1;-1]*semCntrl'),'k','LineWidth',2)
plot(ages(nUcool>=uMin),meanCool(nUcool>=uMin),'c.','MarkerSize',20)
plot(ages(nUcool<uMin),meanCool(nUcool<uMin),'co','MarkerSize',7,'LineWidth',2)
plot(repmat(ages,2,1),meanCool'+([1;-1]*semCool'),'c','LineWidth',2)
xlabel('age (postnatal day)')
ylabel(metric)

if strcmp(metric,'rPref')

    % Plot SI
    figure; hold on
    plot(ages(nUcntrl>=uMin),SI(nUcntrl>=uMin),'k.','MarkerSize',20)
    plot(ages(nUcntrl<uMin),SI(nUcntrl<uMin),'ko','MarkerSize',7,'LineWidth',2)
    plot(repmat(ages,2,1),SI'+([1;-1]*semSI'),'k-','LineWidth',2)
    ageBins = min(ages):3:max(ages)+3;
    for ab = 1:length(ageBins)-1
        ageIdx = ages>=ageBins(ab) & ages<ageBins(ab+1);
        binnedSI(ab) = mean(SI(ageIdx),'omitnan');
        binnedSem(ab) = std(SI(ageIdx),'omitnan')/sqrt(sum(ageIdx));
        binnedAges(ab) = mean(ages(ageIdx),'omitnan');
    end
    nanIdx = isnan(binnedSI);
    binnedSI = binnedSI(~nanIdx);
    binnedSem = binnedSem(~nanIdx);
    binnedAges = binnedAges(~nanIdx);
    p(1) = plot(binnedAges,binnedSI,'r-','LineWidth',2);
    plot(repmat(binnedAges,2,1),binnedSI+([-1;1]*binnedSem),'r-','LineWidth',2)
    [linFit,gofLin] = fit(ages(nID)',SI(nID),'(m*x)+b');
    [poly2Fit,gofPoly2] = fit(ages(nID)',SI(nID),'poly2');
    [poly3Fit,gofPoly3] = fit(ages(nID)',SI(nID),'poly3');
    p(2) = plot(linFit,'r--');
    p(3) = plot(poly2Fit,'g--');
    p(4) = plot(poly3Fit,'b--');
    legend(p,{'binned mean',['linear fit; rmse=' num2str(gofLin.rmse)],['poly 2 fit; rmse=' num2str(gofPoly2.rmse)],['poly 3 fit; rmse=' num2str(gofPoly3.rmse)]})
    xlabel('age (postnatal day)')
    ylabel('SI: (Rcool-Rcntrl)/(Rcool+Rcntrl)')
    
    figure; hold on
    siDist_age = vertcat(siDist_age{:});
    siDist = vertcat(siDist{:});
    plot(siDist_age,siDist,'k.')
    ageBins = min(siDist_age):3:max(siDist_age)+3;
    for ab = 1:length(ageBins)-1
        ageIdx = siDist_age>=ageBins(ab) & siDist_age<ageBins(ab+1);
        binnedSImean(ab) = mean(siDist(ageIdx),'omitnan');
        binnedSIsem(ab) = std(siDist(ageIdx),'omitnan')/sqrt(sum(ageIdx));
        binnedSIage(ab) = mean(unique(siDist_age(ageIdx)));
    end
    nanIdx = isnan(binnedSImean);
    binnedSImean = binnedSImean(~nanIdx);
    binnedSIsem = binnedSIsem(~nanIdx);
    binnedSIage = binnedSIage(~nanIdx);
    plot(binnedSIage,binnedSImean,'r','LineWidth',2)
    plot(repmat(binnedSIage,2,1),binnedSImean+([-1;1]*binnedSIsem),'r','LineWidth',2)
    xlabel('age (postnatal day)')
    ylabel('SI: (Rcool-Rcntrl)/(Rcool+Rcntrl)')

else

    figure; hold on
    plot(ages(nU>=uMin),dMet(nU>=uMin),'k.','MarkerSize',20)
    plot(ages(nU<uMin),dMet(nU<uMin),'ko','MarkerSize',7,'LineWidth',2)
    plot(repmat(ages,2,1),dMet+([1;-1]*sem_dMet),'k','LineWidth',2)
    yline(0,'k--')
    xlabel('age (postnatal day)')
    ylabel(['delta ' metric])
    

end

end





