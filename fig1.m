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
for sh = 1:nShaft
    shaftIdx{sh} = find(sumStats.shaft==sh);
    [~,sortZidx{sh}] = sort(sumStats.zPos(shaftIdx{sh}),'descend');
    shaftIdx{sh} = shaftIdx{sh}(sortZidx{sh});
    shaftIdx{sh} = shaftIdx{sh}(ismember(shaftIdx{sh},find(goodUIdx)));
end

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

% cooling pump information
pumpKey = vertcat(trialId, pumpBit, [1 diff(pumpBit)] == 1, [1 diff(pumpBit)] == -1);
tPumpOff = pumpKey(1,pumpKey(4,:) == 1); tPumpOff = [tPumpOff nTrials];
tPumpOn = pumpKey(1,pumpKey(3,:) == 1); tPumpOn = tPumpOn(2:end);
warmTrials = [];
coldTrials = [];
for t = 1:length(tPumpOff)
    coldTrials = [coldTrials tPumpOff(t)-20:tPumpOff(t)-1];
end
for t = 1:length(tPumpOn)
    warmTrials = [warmTrials tPumpOn(t)-20:tPumpOn(t)-1];
end
blankTrials = sumStats.fr(1).trialNum(:,trialInfo.blankId);
coldTrials = coldTrials(~ismember(coldTrials,blankTrials) & ~(coldTrials<1));
warmTrials = warmTrials(~ismember(warmTrials,blankTrials) & ~(warmTrials<1));

% trial-wise firing rate
for u = 1:nU
    fr(u,:) = sumStats.fr(u).bc(:);
    fr_trialId(u,:) = sumStats.fr(u).trialNum(:);
    [fr_trialId(u,:),sortIdx] = sort(fr_trialId(u,:));
    fr(u,:) = fr(u,sortIdx);
end
fr_norm = fr./max(fr,[],2);

% psth
bw = 6;
edges = (trialStart(1):sf*bw:trialEnd(end))/sf;
for u = 1:nU
    psth(u,:) = histcounts(spks(u).times/sf,edges)/bw;
end
psth_norm = psth./max(psth,[],2);

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
    plot(fr_trialId(1,:),mean(fr(shaftIdx{sh},:),'omitnan'),'o-')
end
ylabel('firing rate (Hz)')
xlim([1 nTrials])
xlabel('trial #')

figure('Position',[100 100 600 1100]); tiledlayout(nShaft+1,1)
for sh = 1:nShaft+1
    nexttile;hold on
    if sh == nShaft+1
        xline(coldTrials,'c--')
        xline(warmTrials,'r--')
        yyaxis left
        plot(trialId,tempDat,'LineWidth',2)
        ylabel('temperature (c)')
        yyaxis right
        plot(trialId,pumpBit,'LineWidth',2)
        ylabel('pump on/off')
        xlabel('trial #')
    else
        imagesc(fr_norm(shaftIdx{sh},:))
        colorbar
        caxis([0 1])
        axis tight
        ylabel(['shaft #' num2str(sh)])
    end
    xlim([1 nTrials])
end

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
    p(sh) = plot(edges(2:end)/60,mean(psth_norm(shaftIdx{sh},:),'omitnan'));
    legLbl{sh} = ['shaft #' num2str(sh)];
end
ylabel('mean normalized response')
xlim([trialStart(1) trialStart(end)]/(sf*60))
xlabel('time (min)')
legend(p,legLbl)

figure('Position',[100 100 600 1100]); tiledlayout(nShaft+1,1);
for sh = 1:nShaft+1
    nexttile; hold on
    if sh == nShaft+1
        xline(trialStart(coldTrials)/(sf*60),'c--')
        xline(trialStart(warmTrials)/(sf*60),'r--')
        yyaxis left
        plot(trialStart(trialId)/(sf*60),tempDat,'LineWidth',2)
        ylabel('temperature (c)')
        yyaxis right
        plot(trialStart(trialId)/(sf*60),pumpBit,'LineWidth',2)
        ylabel('pump on/off')
        xlim([trialStart(1) trialStart(end)]/(sf*60))
        xlabel('time (min)')
    else
        imagesc(psth_norm(shaftIdx{sh},:))
        colorbar
        axis tight
        ylabel(['shaft #' num2str(sh)])
        xt = xticks*(bw/60);
        for tic = 1:length(xt)
            xt_lbl{tic} = num2str(xt(tic));
        end
        xticklabels(xt_lbl)
    end
end

% probe layout
x = sumStats.xPos;
y = sumStats.zPos;
rWarm = mean(fr(:,warmTrials),2,'omitnan');
rCold = mean(fr(:,coldTrials),2,'omitnan');

figure('Position',[100 100 1100 1100]); tiledlayout(1,2)
rLims = [0 70];
nexttile; hold on
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


