clear all
close all

animal = 'febq5';
unit = '001';
expt = '007';
probe = 1;
anaMode = 'MU';
% dataFold = 'F:\Brandon\data';
dataFold = 'Y:\Brandon\data';
% dataFold = 'C:\Users\brand\Documents\data';
% dataFold = '/Volumes/Lab drive/Brandon/data';
% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
exptName = [animal '_u' unit '_' expt];

% load necessary data
load(fullfile(dataFold,'Ephys',animal,exptName,[exptName '_temp.mat']))
load(fullfile(dataFold,'Ephys',animal,exptName,[exptName '_trialInfo.mat']))
load(fullfile(dataFold,'Ephys',animal,exptName,[exptName '_id.mat']))
load(fullfile(dataFold,'Ephys',animal,exptName,[exptName '.analyzer']),'-mat')
[sumStats,spks] = anaOri(animal,unit,expt,probe,anaMode,dataFold,0,0);
nU = height(sumStats);
nShaft = length(unique(sumStats.shaft));
goodUIdx = screenUnits(sumStats,'MU');

% Organize groups by shaft
shaftIdx = sumStats.shaft;

% Organize groups by depth
depthLims = min(sumStats.zPos):100:max(sumStats.zPos);
depthIdx = zeros(1,nU);
for d = 1:(length(depthLims)-1)
    depthIdx(sumStats.zPos>=depthLims(d) & sumStats.zPos<depthLims(d+1)) = d;
end
depthIdx(sumStats.zPos>=depthLims(d+1) & sumStats.zPos<=max(sumStats.zPos)) = d+1;

% trial information
pre = getparam('predelay',Analyzer);
stimTime = getparam('stim_time',Analyzer);
post = getparam('postdelay',Analyzer);

eventIDdiff = [1;diff(trialInfo.eventId)];
eventTimes = trialInfo.eventTimes;
trialStart = eventTimes(eventIDdiff==1);
stimStart = eventTimes(eventIDdiff==2);
stimEnd = eventTimes(eventIDdiff==-2);
trialEnd = eventTimes(eventIDdiff==-1);

sf = id.sampleFreq;
nTrials = length(trialInfo.triallist);
blankTrials = sumStats.fr(1).trialNum(:,trialInfo.blankId);

% cooling pump information
warmTrials = [];
warmT = trialId(tempDat>30);
gapInd = [0 diff(warmT)]>5;
warmStart = [warmT(1) warmT(gapInd)];
warmEnd = [warmT(find(gapInd)-1) warmT(end)];
for t = 1:length(warmStart)
    warmTrials = [warmTrials warmStart(t):warmEnd(t)];
end
coldTrials = [];
coldT = trialId(tempDat<8);
gapInd = [0 diff(coldT)]>5;
coldStart = [coldT(1) coldT(gapInd)];
coldEnd = [coldT(find(gapInd)-1) coldT(end)];
for t = 1:length(coldStart)
    coldTrials = [coldTrials coldStart(t):coldEnd(t)];
end
pumpKey = vertcat(trialId, pumpBit, [1 diff(pumpBit)] == 1, [1 diff(pumpBit)] == -1);
tPumpOff = pumpKey(1,pumpKey(4,:) == 1); tPumpOff = [tPumpOff nTrials];
tPumpOn = pumpKey(1,pumpKey(3,:) == 1); tPumpOn = tPumpOn(2:end);
pumpSqWv = [0;1];
for i = 1:length(tPumpOn)
    pumpSqWv = [pumpSqWv,[0;tPumpOn(i)],[1;tPumpOn(i)],[1;tPumpOff(i)],[0;tPumpOff(i)]];
end
% warmTrials = [];
% coldTrials = [];
% for t = 1:length(tPumpOff)
%     coldTrials = [coldTrials tPumpOff(t)-32:tPumpOff(t)-1];
% end
% for t = 1:length(tPumpOn)
%     warmTrials = [warmTrials tPumpOn(t)-32:tPumpOn(t)-1];
% end
coldTrials = coldTrials(~ismember(coldTrials,blankTrials) & ~(coldTrials<1));
warmTrials = warmTrials(~ismember(warmTrials,blankTrials) & ~(warmTrials<1));

% trial-wise firing rate
for u = 1:nU
    prefC = find(sumStats.condition{u}==sumStats.oriPref(u));
    prefIdx = sumStats.fr(u).trialCond==prefC;
    fr_pref(u,:) = sumStats.fr(u).stim(prefIdx);
    fr_pref_trialId(u,:) = sumStats.fr(u).trialNum(prefIdx);
    fr(u,:) = sumStats.fr(u).stim(:);
    fr_trialId(u,:) = sumStats.fr(u).trialNum(:);
    [fr_trialId(u,:),sortIdx] = sort(fr_trialId(u,:));
    fr(u,:) = fr(u,sortIdx);
end
fr_norm = fr./max(fr,[],2);

% psth
bw = 12;
edges = (trialStart(1):sf*bw:trialEnd(end))/sf;
for u = 1:nU
    psth(u,:) = histcounts(spks(u).times/sf,edges)/bw;
end
psth_norm = psth./max(psth,[],2);

for u = 1:nU
    warmTrialIdx = ismember(fr_pref_trialId(u,:),warmTrials);
    rWarm(u) = mean(fr_pref(u,warmTrialIdx),'omitnan');
    nWarmTrials(u) = sum(warmTrialIdx);
    coldTrialIdx = ismember(fr_pref_trialId(u,:),coldTrials);
    rCold(u) = mean(fr_pref(u,coldTrialIdx),'omitnan');
    nColdTrials(u) = sum(coldTrialIdx);
end
si = (rCold-rWarm)./(rCold+rWarm);

%% plot

% stim period (baseline corrected) firing rate
figure; hold on
yyaxis left
plot(fr_trialId(1,:),mean(fr(goodUIdx,:),'omitnan'),'o-')
ylim([0 18])
ylabel('firing rate (Hz)')
yyaxis right
plot(trialId,tempDat,'LineWidth',2)
ylabel('temperature (C)')
xlim([1 nTrials])
xlabel('trial #')

figure; hold on
for sh = 1:nShaft
    plot(fr_trialId(1,:),mean(fr(shaftIdx==sh,:),'omitnan'),'o-')
end
ylabel('firing rate (Hz)')
xlim([1 nTrials])
xlabel('trial #')

figure; tiledlayout(nShaft+2,1)
nexttile; hold on
xline(coldTrials,'c--')
xline(warmTrials,'r--')
plot(trialId,tempDat,'k','LineWidth',2)
ylabel('temperature (c)')
xlim([1 nTrials])
nexttile; hold on
%         plot(trialId,pumpBit,'LineWidth',2)
plot(pumpSqWv(2,:),pumpSqWv(1,:),'k','LineWidth',2)
ylabel('pump on/off')
xlabel('trial number')
xlim([1 nTrials])

% figure; tiledlayout(nShaft,1)
for sh = 1:nShaft
    nexttile;hold on
    [~,zPosSort] = sort(sumStats.zPos(shaftIdx==sh),'descend');
    r = fr_norm(shaftIdx==sh,:);
    imagesc(r(zPosSort,:))
    colormap gray
    colorbar
    caxis([0 1])
    axis tight
    ylabel(['shaft #' num2str(sh)])
    xlim([1 nTrials])
end
xlabel('trial number')

% psth
figure; hold on
yyaxis left
plot(edges(2:end)/60,mean(psth_norm(goodUIdx,:),'omitnan'))
ylabel('mean normalized response')
yyaxis right
plot(trialStart(trialId)/(sf*60),tempDat,'LineWidth',2)
ylabel('temperature (C)')
xlim([trialStart(1) trialStart(end)]/(sf*60))
xlabel('time (min)')

figure; hold on
for sh = 1:nShaft
    p(sh) = plot(edges(2:end)/60,mean(psth_norm(shaftIdx==sh,:),'omitnan'));
    legLbl{sh} = ['shaft #' num2str(sh)];
end
ylabel('mean normalized response')
xlim([trialStart(1) trialStart(end)]/(sf*60))
xlabel('time (min)')
legend(p,legLbl)

figure; tiledlayout(nShaft+1,1);
for sh = 1:nShaft+1
    nexttile; hold on
    if sh == nShaft+1
        xline(trialStart(coldTrials)/(sf*60),'c--')
        xline(trialStart(warmTrials)/(sf*60),'r--')
        yyaxis left
        plot(trialStart(trialId)/(sf*60),tempDat,'LineWidth',2)
        ylabel('temperature (c)')
        yyaxis right
%         plot(trialStart(trialId)/(sf*60),pumpBit,'LineWidth',2)
        plot(trialStart(pumpSqWv(2,:))/(sf*60),pumpSqWv(1,:),'LineWidth',2)
        ylabel('pump on/off')
        xlim([trialStart(1) trialStart(end)]/(sf*60))
        xlabel('time (min)')
    else
        [~,zPosSort] = sort(sumStats.zPos(shaftIdx==sh),'descend');
        r = psth_norm(shaftIdx==sh,:);
        imagesc(r(zPosSort,:))
        colormap gray
        colorbar
        axis tight
        ylabel(['shaft #' num2str(sh)])
        xticks((5:5:30)/(bw/60))
        xt = xticks*(bw/60);
        for tic = 1:length(xt)
            xt_lbl{tic} = num2str(xt(tic));
        end
        xticklabels(xt_lbl)
    end
end

% probe layout

figure
tiledlayout(1,2)
rLims = [0 70];
nexttile; hold on
x = sumStats.xPos;
y = sumStats.zPos;
bubblechart(x, y, rWarm, rWarm)
bubblelim(rLims)
caxis(rLims)
scatter(x(~goodUIdx), y(~goodUIdx),'rx')
set(gca,'YDir','reverse')
nexttile; hold on
bubblechart(x, y, rCold, rCold)
bubblelegend('firing rate','Location','eastoutside')
bubblelim(rLims)
caxis(rLims)
colorbar
scatter(x(~goodUIdx), y(~goodUIdx),'rx')
set(gca,'YDir','reverse')

figure;
chs = 97;%(depthIdx==1|depthIdx==2) & shaftIdx'==4;
subplot(1,2,2);hold on
x = sumStats.xPos;
y = sumStats.zPos;
z = si;
scatter(x, y,'k')
scatter(x, y,[],z,'filled')
scatter(x(chs),y(chs),'g')
caxis([-1 1])
colormap cool
cb = colorbar;
cb.Label.String = 'SI';
cb.Label.FontSize = 10;
% scatter(x(~goodUIdx), y(~goodUIdx),'rx')
set(gca,'YDir','reverse')
xlim([-200 1200])
ylabel('depth')
xlabel('distance from cooling loop')
box on
subplot(2,2,1);hold on
spks = [sumStats.spkTimes{chs}];
patch([0 1 1 0],[0 0 nTrials+1 nTrials+1],'k','EdgeColor','none','FaceAlpha',0.2)
for i = 1:length(coldStart)
patch([-1 2 2 -1],[coldStart(i) coldStart(i) coldEnd(i) coldEnd(i)],'c','EdgeColor','none','FaceAlpha',0.2)
end
for i = 1:length(warmStart)
patch([-1 2 2 -1],[warmStart(i) warmStart(i) warmEnd(i) warmEnd(i)],'r','EdgeColor','none','FaceAlpha',0.2)
end
plot(spks(1,:),spks(2,:),'k.','MarkerSize',3)
ylim([0 nTrials+1])
ylabel('trial number')
xlim([-1 2])
xlabel('time (sec)')
box on
axis square
subplot(2,2,3);hold on
x = edges(2:end)/60;
y = psth(chs,:);
plot(x,y,'k')
ylabel('mean norm. psth')
yyaxis right
plot(trialStart(trialId)/(sf*60),tempDat,'LineWidth',2)
ylabel('temperature (C)')
xlabel('time (min)')
% xlim([0 nTrials])
box on
axis square

figure; hold on
cntrlStats = anaOri('febq5','001','004',1,'MU',dataFold,0,0);
x = cntrlStats.rPref;
axLim = [0 80];
scatter(x,rCold,'c.')
scatter(x(chs),rCold(chs),'go')
fit1 = polyfit(x,rCold,1);
plot(axLim,polyval(fit1,axLim),'c-')
scatter(x,rWarm,'k.')
scatter(x(chs),rWarm(chs),'go')
fit2 = polyfit(x,rWarm,1);
plot(axLim,polyval(fit2,axLim),'k-')
plot(axLim,axLim,'k--')
xlim(axLim)
xlabel('pre-cooling firing rate (Hz)')
ylim(axLim)
ylabel('firing rate (Hz)')
box on
axis square


% figure; hold on
% cntrlStats = anaOri('febq5','001','004',1,'MU',dataFold,0,0);
% x = cntrlStats.rPref;
% axLim = [0 80];
% shapes = {'o','square','diamond','^'};
% for sh = unique(shaftIdx)'
%     idx = shaftIdx == sh;
%     scatter(x(idx),rCold(idx),['c' shapes{sh}])
%     fit1 = polyfit(x(idx),rCold(idx),1);
%     plot(axLim,polyval(fit1,axLim),'c-')
%     scatter(x(idx),rWarm(idx),['k' shapes{sh}])
%     fit2 = polyfit(x(idx),rWarm(idx),1);
%     plot(axLim,polyval(fit2,axLim),'k-')
% end
% plot(axLim,axLim,'k--')
% xlim(axLim)
% ylim(axLim)
% box on
% axis square





% % for febq5_u001_007: plot disinhibited channels on shaft 1
% figure;
% subplot(1,2,1); hold on
% scatter(sumStats.xPos, sumStats.zPos,'k')
% scatter(sumStats.xPos, sumStats.zPos,[],si,'filled')
% patch([900 1100 1100 900],[600 600 400 400],'r','EdgeColor','r','FaceColor','none')
% colormap cool
% cb = colorbar;
% cb.Label.String = 'SI';
% cb.Label.FontSize = 10;
% set(gca,'YDir','reverse')
% xlim([-200 1200])
% ylabel('depth')
% xlabel('distance from cooling loop')
% box on
% subplot(2,2,2); hold on
% chs = find(sumStats.xPos>800 & (sumStats.zPos>400 & sumStats.zPos<600));
% spks = [sumStats.spkTimes{chs}];
% for i = 1:length(coldStart)
% patch([-1 2 2 -1],[coldStart(i) coldStart(i) coldEnd(i) coldEnd(i)],'c','EdgeColor','none','FaceAlpha',0.2)
% end
% for i = 1:length(warmStart)
% patch([-1 2 2 -1],[warmStart(i) warmStart(i) warmEnd(i) warmEnd(i)],'r','EdgeColor','none','FaceAlpha',0.2)
% end
% plot(spks(1,:),spks(2,:),'k.','MarkerSize',3)
% ylim([0 nTrials+1])
% ylabel('trial number')
% xlim([-1 2])
% xlabel('time (sec)')
% box on
% subplot(2,2,4); hold on
% y = mean(psth_norm(chs,:),'omitnan');
% x = edges(2:end)/60;
% plot(x,y,'k')
% ylabel('mean norm. psth')
% yyaxis right
% plot(trialStart(trialId)/(sf*60),tempDat,'LineWidth',2)
% ylabel('temperature (C)')
% xlabel('time (min)')
% box on
% sgtitle('febq5 001 007: disinhibited channels')










% for u = goodUnitId'
% 
% figure;
% tiledlayout(2,2)
% bw = id.sampleFreq*1;
% 
% uIdx = [spks.unitId] == u;
% nexttile; hold on
% spkTimes = spks(uIdx).stimCent;
% patch([0 1 1 0],[0 0 nTrials+1 nTrials+1],'k','EdgeColor','none','FaceAlpha',0.2)
% scatter(spkTimes(1,:),spkTimes(2,:),1,'k.')
% ylim([0 nTrials+1])
% ylabel('trial')
% xlim([-1 2])
% xlabel('time (sec)')
% title(['channel ' num2str(spks(uIdx).unitId)])
% box on
% 
% nexttile; hold on
% colororder({'k','r'})
% [psth,edges] = histcounts(spks(uIdx).times,'BinWidth',bw);
% yyaxis left
% plot(edges(2:end)/(sf*60),psth)
% ylabel('firing rate (hz)')
% yyaxis right
% plot(trialStart(trialId)/(sf*60),tempDat)
% ylabel('temperature (c)')
% xlabel('time (min)')
% box on
% 
% nexttile; hold on
% spkTimes = [spks(:).stimCent];
% patch([0 1 1 0],[0 0 nTrials+1 nTrials+1],'k','EdgeColor','none','FaceAlpha',0.2)
% scatter(spkTimes(1,:),spkTimes(2,:),1,'k.')
% ylim([0 nTrials+1])
% ylabel('trial')
% xlim([-1 2])
% xlabel('time (sec)')
% title("'good' channels")
% box on
% 
% nexttile; hold on
% [psth,edges] = histcounts([spks.times],'BinWidth',bw);
% psth = psth;
% yyaxis left
% plot(edges(2:end)/(sf*60),psth)
% ylabel('firing rate (hz)')
% yyaxis right
% plot(trialStart(trialId)/(sf*60),tempDat)
% ylabel('')
% xlabel('time (min)')
% box on
% 
% sgtitle([animal ' ' unit ' ' expt])
% 
% end
% 
% % figure; hold on
% % scatter([spks.xPos],[spks.zPos])
% % scatter([spks(~goodIdx).xPos],[spks(~goodIdx).zPos],'rx')
% % scatter(sumStats.xPos(goodIdx),sumStats.zPos(goodIdx),200,repmat( warmR', 1,3 ),'.')
% % 
% % figure; hold on
% % scatter([spks.xPos],[spks.zPos])
% % scatter([spks(~goodIdx).xPos],[spks(~goodIdx).zPos],'rx')
% % scatter(sumStats.xPos(goodIdx),sumStats.zPos(goodIdx),200,repmat( coldR', 1,3 ),'.')


