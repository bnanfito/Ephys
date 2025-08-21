%anaV1cool_ori2
clear
close all

anaMode = 'MU';
proj = 'V1cool_ori';
area = 'PSS';
dataFold = 'Y:\Brandon\data';
% dataFold = '/Volumes/NielsenHome2/Brandon/data';
% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
ageGroups = {[28 32],[33 40],[41 150]};

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
for ag = 1:length(ageGroups)
    curAGidx = find(ages>=ageGroups{ag}(1) & ages<=ageGroups{ag}(2));
    agIdx(curAGidx) = ag;
end
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

metric = 'latency';

SI = nan(length(animals),1);
semSI = nan(length(animals),1);
meanCntrl = nan(length(animals),1);
semCntrl = nan(length(animals),1);
meanCool = nan(length(animals),1);
semCool = nan(length(animals),1);
for a = 1:length(animals)
    if ~(isempty(dat.cntrl{a,1}) || isempty(dat.cool{a,1}))

        rPref.cntrl.dist{a} = dat.cntrl{a,1}.rPref(goodIdCntrl{a,1});
        rPref.cntrl.ave(a) = mean(rPref.cntrl.dist{a},'omitnan');
        rPref.cntrl.sem(a) = std(rPref.cntrl.dist{a},'omitnan')/sqrt(length(rPref.cntrl.dist{a}));
        rPref.cntrl.n(a) = length(rPref.cntrl.dist{a});

        rPref.cool.dist{a} = dat.cool{a,1}.rPref(goodIdCntrl{a,1});
        rPref.cool.ave(a) = mean(rPref.cool.dist{a},'omitnan');
        rPref.cool.sem(a) = std(rPref.cool.dist{a},'omitnan')/sqrt(length(rPref.cool.dist{a}));
        rPref.cool.n(a) = length(rPref.cool.dist{a});

        si.dist{a} = (rPref.cool.dist{a}-rPref.cntrl.dist{a})./(rPref.cool.dist{a}+rPref.cntrl.dist{a});
        si.ave(a) = mean(si.dist{a},'omitnan');
        si.sem(a) = std(si.dist{a},'omitnan')/sqrt(length(si.dist{a}));
        si.n(a) = length(si.dist{a});

        late.cntrl.dist{a} = dat.cntrl{a,1}.latency(goodIdCntrl{a,1});
        late.cntrl.ave(a) = mean(late.cntrl.dist{a},'omitnan');
        late.cntrl.sem(a) = std(late.cntrl.dist{a},'omitnan')/sqrt(length(late.cntrl.dist{a}));
        late.cntrl.n(a) = length(late.cntrl.dist{a});

        late.cool.dist{a} = dat.cool{a,1}.latency(goodIdCntrl{a,1});
        late.cool.ave(a) = mean(late.cool.dist{a},'omitnan');
        late.cool.sem(a) = std(late.cool.dist{a},'omitnan')/sqrt(length(late.cool.dist{a}));
        late.cool.n(a) = length(late.cool.dist{a});

        dLate.dist{a} = (late.cool.dist{a}-late.cntrl.dist{a})./(late.cool.dist{a}+late.cntrl.dist{a});
        dLate.ave(a) = mean(dLate.dist{a},'omitnan');
        dLate.sem(a) = std(dLate.dist{a},'omitnan')/sqrt(length(dLate.dist{a}));
        dLate.n(a) = length(dLate.dist{a});

%         % Distributions
%         if strcmp(metric,'response var')
%             cntrlDist = squeeze(mean(std(cat(3,dat.cntrl{a,1}.response{:}),'omitnan'),'omitnan'));
%             coolDist = squeeze(mean(std(cat(3,dat.cool{a,1}.response{:}),'omitnan'),'omitnan'));
%         else
%             cntrlDist = dat.cntrl{a,1}{:,strcmp(dat.cntrl{a,1}.Properties.VariableNames,metric)}(goodIdCntrl{a,1});
%             coolDist = dat.cool{a,1}{:,strcmp(dat.cool{a,1}.Properties.VariableNames,metric)}(goodIdCntrl{a,1});
%         end
%         meanCntrl(a) = mean(cntrlDist,'omitnan');
%         semCntrl(a) = std(cntrlDist,'omitnan')/sqrt(length(cntrlDist));
%         nUcntrl(a,1) = length(cntrlDist);
%         meanCool(a) = mean(coolDist,'omitnan');
%         semCool(a) = std(coolDist,'omitnan')/sqrt(length(coolDist));
%         nUcool(a,1) = length(coolDist);
% 
%         if strcmp(metric,'rPref')
%             % Suppression Index
%             siDist{a} = (coolDist-cntrlDist)./(coolDist+cntrlDist);
%             siDist_age{a} = repmat(ages(a),length(siDist{a}),1);
%             SI(a) = mean(siDist{a},'omitnan');
%             semSI(a) = std(siDist{a},'omitnan')/sqrt(length(siDist{a}));
%         else
%             % Delta Ldir
%             goodId_cntrl = dat.cntrl{a,1}.uID(screenUnits(dat.cntrl{a,1},anaMode));
%             goodId_cool = dat.cool{a,1}.uID(screenUnits(dat.cool{a,1},anaMode));
%             goodId = goodId_cntrl(ismember(goodId_cntrl,goodId_cool));
%             nU(a,1) = length(goodId);
%             dMetDist = coolDist-cntrlDist;
%             dMet(a) = mean(dMetDist,'omitnan');
%             sem_dMet(a) = std(dMetDist,'omitnan')/sqrt(length(dMetDist));
%         end

    else
        rPref.cntrl.dist{a} = [];
        rPref.cntrl.ave(a) = nan;
        rPref.cntrl.sem(a) = nan;
        rPref.cntrl.n(a) = nan;

        rPref.cool.dist{a} = [];
        rPref.cool.ave(a) = nan;
        rPref.cool.sem(a) = nan;
        rPref.cool.n(a) = nan;

        si.dist{a} = [];
        si.ave(a) = nan;
        si.sem(a) = nan;
        si.n(a) = nan;

        late.cntrl.dist{a} = [];
        late.cntrl.ave(a) = nan;
        late.cntrl.sem(a) = nan;
        late.cntrl.n(a) = nan;

        late.cool.dist{a} = [];
        late.cool.ave(a) = nan;
        late.cool.sem(a) = nan;
        late.cool.n(a) = nan;

        dLate.dist{a} = [];
        dLate.ave(a) = nan;
        dLate.sem(a) = nan;
        dLate.n(a) = nan;
    end
end
nID = ~isnan(SI);

%% Plot

% a = 4; 
% for u = [42 52 54]
% d1 = dat.cntrl{a}.spkTimes{u}(1,:);
% d2 = dat.cool{a}.spkTimes{u}(1,:);
% nTrials1 = max(dat.cntrl{a}.fr(u).trialNum,[],'all');
% nTrials2 = max(dat.cool{a}.fr(u).trialNum,[],'all');
% bw = 0.01;
% bins = -1:bw:2;
% h1 = histcounts(d1,bins)/(nTrials1*bw);
% h2 = histcounts(d2,bins)/(nTrials2*bw);
% figure; hold on
% maxH = 60;%max([h1 h2],[],'all')+1;
% patch([0 1 1 0],[0 0 maxH maxH],'k','EdgeColor','none','FaceAlpha',0.2)
% plot(bins(2:end),h1,'k')
% plot(bins(2:end),h2,'c')
% xlim([-1 2])
% ylim([0 maxH])
% box on
% axis square
% end


figure
subplot(2,2,1); hold on
idx = agIdx==1|agIdx==2;
xLims = [min(ages(idx))-2 max(ages(idx))+2 max(ages(idx))+2 min(ages(idx))-2];
ag3 = rPref.cntrl.ave(agIdx==3);
mAg3 = mean(ag3,'omitnan');
semAg3 = std(ag3,'omitnan')/sqrt(length(ag3));
y1 = mAg3-semAg3;
y2 = mAg3+semAg3;
y = [y1 y1 y2 y2];
patch(xLims,y,'k','EdgeColor','none','FaceAlpha',0.2)
yline(mAg3,'k--')
ag3 = rPref.cool.ave(agIdx==3);
mAg3 = mean(ag3,'omitnan');
semAg3 = std(ag3,'omitnan')/sqrt(length(ag3));
y1 = mAg3-semAg3;
y2 = mAg3+semAg3;
y = [y1 y1 y2 y2];
patch(xLims,y,'c','EdgeColor','none','FaceAlpha',0.2)
yline(mAg3,'c--')
errorbar(ages(idx),rPref.cntrl.ave(idx),rPref.cntrl.sem(idx),'k','LineStyle','none')
errorbar(ages(idx),rPref.cool.ave(idx),rPref.cool.sem(idx),'c','LineStyle','none')
xlabel('age')
xlim(xLims(1:2))
ylabel('firing rate')
box on

subplot(2,2,2); hold on
ag3 = si.ave(agIdx==3);
mAg3 = mean(ag3,'omitnan');
semAg3 = std(ag3,'omitnan')/sqrt(length(ag3));
y1 = mAg3-semAg3;
y2 = mAg3+semAg3;
y = [y1 y1 y2 y2];
patch(xLims,y,'r','EdgeColor','none','FaceAlpha',0.2)
yline(mAg3,'r--')
yline(0,'k--')
errorbar(ages(idx),si.ave(idx),si.sem(idx),'LineStyle','none')
xlabel('age')
xlim(xLims(1:2))
ylabel('si')
box on

subplot(2,2,3); hold on
ag3 = late.cntrl.ave(agIdx==3);
mAg3 = mean(ag3,'omitnan');
semAg3 = std(ag3,'omitnan')/sqrt(length(ag3));
y1 = mAg3-semAg3;
y2 = mAg3+semAg3;
y = [y1 y1 y2 y2];
patch(xLims,y,'k','EdgeColor','none','FaceAlpha',0.2)
yline(mAg3,'k--')
ag3 = late.cool.ave(agIdx==3);
mAg3 = mean(ag3,'omitnan');
semAg3 = std(ag3,'omitnan')/sqrt(length(ag3));
y1 = mAg3-semAg3;
y2 = mAg3+semAg3;
y = [y1 y1 y2 y2];
patch(xLims,y,'c','EdgeColor','none','FaceAlpha',0.2)
yline(mAg3,'c--')
errorbar(ages(idx),late.cntrl.ave(idx),late.cntrl.sem(idx),'k','LineStyle','none')
errorbar(ages(idx),late.cool.ave(idx),late.cool.sem(idx),'c','LineStyle','none')
xlabel('age')
xlim(xLims(1:2))
ylabel('latency')
box on

subplot(2,2,4); hold on
ag3 = dLate.ave(agIdx==3);
mAg3 = mean(ag3,'omitnan');
semAg3 = std(ag3,'omitnan')/sqrt(length(ag3));
y1 = mAg3-semAg3;
y2 = mAg3+semAg3;
y = [y1 y1 y2 y2];
patch(xLims,y,'r','EdgeColor','none','FaceAlpha',0.2)
yline(mAg3,'r--')
yline(0,'k--')
errorbar(ages(idx),dLate.ave(idx),dLate.sem(idx),'LineStyle','none')
xlabel('age')
xlim(xLims(1:2))
ylabel('delta latency')
box on







% uMin = 20;
% scl = 10;
% 
% if strcmp(anaMode,'MU')
% 
% % Plot distributions
% figure('Position',[100 100 600 600]); hold on
% idx1 = nUcntrl>=uMin & (agIdx(:,1) | agIdx(:,2));
% idx2 = nUcntrl<uMin & (agIdx(:,1) | agIdx(:,2));
% mAG3_cntrl = mean(meanCntrl(agIdx(:,3)),'omitnan');
% semAG3_cntrl = std(meanCntrl(agIdx(:,3)),'omitnan')/sqrt(sum(agIdx(:,3)));
% mAG3_cool = mean(meanCool(agIdx(:,3)),'omitnan');
% semAG3_cool = std(meanCool(agIdx(:,3)),'omitnan')/sqrt(sum(agIdx(:,3)));
%     x = [min(ages(idx1|idx2))-1  max(ages(idx1|idx2))+1 max(ages(idx1|idx2))+1  min(ages(idx1|idx2))-1];
%     y = [mAG3_cntrl-semAG3_cntrl mAG3_cntrl-semAG3_cntrl mAG3_cntrl+semAG3_cntrl mAG3_cntrl+semAG3_cntrl];
%     patch(x,y,'k','EdgeColor','none','FaceAlpha',0.2)
%     yline(mAG3_cntrl,'k--','LineWidth',2)
%     y = [mAG3_cool-semAG3_cool mAG3_cool-semAG3_cool mAG3_cool+semAG3_cool mAG3_cool+semAG3_cool];
%     patch(x,y,'c','EdgeColor','none','FaceAlpha',0.2)
%     yline(mAG3_cool,'c--','LineWidth',2)
% errorbar(ages(~isnan(meanCntrl)), meanCntrl(~isnan(meanCntrl)), semCntrl(~isnan(meanCntrl)),'k','LineWidth',2, 'LineStyle','none')
% scatter(ages(idx1 & ~isnan(meanCntrl)),meanCntrl(idx1 & ~isnan(meanCntrl)),10*scl,'k','filled','o','MarkerFaceAlpha',1,'MarkerEdgeColor','k')
% scatter(ages(idx2 & ~isnan(meanCntrl)),meanCntrl(idx2 & ~isnan(meanCntrl)),10*scl,'k','o','LineWidth',2)
% errorbar(ages(~isnan(meanCool)), meanCool(~isnan(meanCool)), semCool(~isnan(meanCool)),'c','LineWidth',2, 'LineStyle','none')
% scatter(ages(idx1 & ~isnan(meanCool)),meanCool(idx1 & ~isnan(meanCool)),10*scl,'c','filled','o','MarkerFaceAlpha',1,'MarkerEdgeColor','k')
% scatter(ages(idx2 & ~isnan(meanCool)),meanCool(idx2 & ~isnan(meanCool)),10*scl,'c','o','LineWidth',2)
% xlabel('age (postnatal day)')
% xlim([min(x) max(x)])
% ylabel(metric)
% box on
% 
% if strcmp(metric,'rPref')
% 
%     % Plot SI
%     figure('Position',[100 100 600 600]); hold on
%     mAG3 = mean(SI(agIdx(:,3)),'omitnan');
%     semAG3 = std(SI(agIdx(:,3)),'omitnan')/sqrt(sum(agIdx(:,3)));
%     x = [min(ages(idx1|idx2))-1  max(ages(idx1|idx2))+1 max(ages(idx1|idx2))+1  min(ages(idx1|idx2))-1];
%     y = [mAG3-semAG3 mAG3-semAG3 mAG3+semAG3 mAG3+semAG3];
%     patch(x,y,'r','EdgeColor','none','FaceAlpha',0.2)
%     yline(mAG3,'r--','LineWidth',2)
%     yline(0,'k--','LineWidth',2)
%     errorbar(ages(~isnan(SI)), SI(~isnan(SI)), semSI(~isnan(SI)),'k','LineWidth',2,'LineStyle','none')
%     scatter(ages(idx1 & ~isnan(SI)),SI(idx1 & ~isnan(SI)),10*scl,'k','filled','o','MarkerFaceAlpha',1,'MarkerEdgeColor','k')
%     scatter(ages(idx2 & ~isnan(SI)),SI(idx2 & ~isnan(SI)),10*scl,'k','o','LineWidth',2)
%     xlabel('age (postnatal day)')
%     xlim([min(x) max(x)])
%     ylabel('SI: (Rcool-Rcntrl)/(Rcool+Rcntrl)')
%     box on
% 
%     
% %     figure; hold on
% %     siDist_age = vertcat(siDist_age{:});
% %     siDist = vertcat(siDist{:});
% %     plot(siDist_age,siDist,'k.')
% %     ageBins = min(siDist_age):3:max(siDist_age)+3;
% %     for ab = 1:length(ageBins)-1
% %         ageIdx = siDist_age>=ageBins(ab) & siDist_age<ageBins(ab+1);
% %         binnedSImean(ab) = mean(siDist(ageIdx),'omitnan');
% %         binnedSIsem(ab) = std(siDist(ageIdx),'omitnan')/sqrt(sum(ageIdx));
% %         binnedSIage(ab) = mean(unique(siDist_age(ageIdx)));
% %     end
% %     nanIdx = isnan(binnedSImean);
% %     binnedSImean = binnedSImean(~nanIdx);
% %     binnedSIsem = binnedSIsem(~nanIdx);
% %     binnedSIage = binnedSIage(~nanIdx);
% %     plot(binnedSIage,binnedSImean,'r','LineWidth',2)
% %     plot(repmat(binnedSIage,2,1),binnedSImean+([-1;1]*binnedSIsem),'r','LineWidth',2)
% %     xlabel('age (postnatal day)')
% %     ylabel('SI: (Rcool-Rcntrl)/(Rcool+Rcntrl)')
% 
% else
% 
%     figure; hold on
%     plot(ages(nU>=uMin),dMet(nU>=uMin),'k.','MarkerSize',20)
%     plot(ages(nU<uMin),dMet(nU<uMin),'ko','MarkerSize',7,'LineWidth',2)
%     plot(repmat(ages,2,1),dMet+([1;-1]*sem_dMet),'k','LineWidth',2)
%     yline(0,'k--')
%     xlabel('age (postnatal day)')
%     ylabel(['delta ' metric])
%     
% 
% end
% 
% end





