clear all
close all

animal = 'febq5';
unit = '001';
expt = '007';
probe = 1;
anaMode = 'MU';
% dataFold = 'F:\Brandon\data';
% dataFold = 'C:\Users\brand\Documents\data';
% dataFold = '/Volumes/Lab drive/Brandon/data';
dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
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
tPumpOff = pumpKey(1,pumpKey(4,:) == 1);
tPumpOn = pumpKey(1,pumpKey(3,:) == 1); tPumpOn = tPumpOn(2:end);
warmTrials = [];
coldTrials = [];
for t = 1:length(tPumpOff)
    coldTrials = [coldTrials tPumpOff(t)-16:tPumpOff(t)-1];
end
for t = 1:length(tPumpOn)
    warmTrials = [warmTrials tPumpOn(t)-16:tPumpOn(t)-1];
end
blankTrials = sumStats.fr(1).trialNum(:,trialInfo.blankId);
coldTrials = coldTrials(~ismember(coldTrials,blankTrials));
warmTrials = warmTrials(~ismember(warmTrials,blankTrials));

% trial-wise firing rate
for u = 1:nU
    fr(u,:) = sumStats.fr(u).bc(:);
    fr_trialId(u,:) = sumStats.fr(u).trialNum(:);
    [fr_trialId(u,:),sortIdx] = sort(fr_trialId(u,:));
    fr(u,:) = fr(u,sortIdx);
end
fr_norm = fr./max(fr,[],2);

% psth
bw = 30;
edges = (trialStart(1):sf*bw:trialEnd(end))/sf;
for u = 1:nU
    psth(u,:) = histcounts(spks(u).times/sf,edges)/(nTrials*bw);
end
psth_norm = psth./max(psth,[],2);

%% plot

figure; hold on
plot(fr_trialId(1,:),mean(fr(goodUIdx,:),'omitnan'),'o-')

figure; hold on
for sh = 1:nShaft
    plot(fr_trialId(1,:),mean(fr(shaftIdx{sh},:),'omitnan'),'o-')
end

figure; tiledlayout(nShaft+1,1)
for sh = 1:nShaft+1
    nexttile;hold on
    if sh == nShaft+1
        yyaxis left
        plot(trialId,tempDat)
        ylabel('temperature (c)')
        yyaxis right
        plot(trialId,pumpBit)
        ylabel('pump on/off')
        xlabel('trial #')
    else
        imagesc(fr_norm(shaftIdx{sh},:))
        colorbar
        axis tight
        ylabel(['shaft #' num2str(sh)])
    end
    xlim([1 nTrials])
end

figure; hold on
plot(edges(2:end),mean(psth_norm(goodUIdx,:),'omitnan'))

figure; hold on
for sh = 1:nShaft
    p(sh) = plot(edges(2:end),mean(psth_norm(shaftIdx{sh},:),'omitnan'));
    legLbl{sh} = ['shaft #' num2str(sh)];
end
legend(p,legLbl)

figure; tiledlayout(nShaft+1,1);
for sh = 1:nShaft+1
    nexttile; hold on
    if sh == nShaft+1
        yyaxis left
        plot(trialStart(trialId)/(sf*60),tempDat)
        ylabel('temperature (c)')
        yyaxis right
        plot(trialStart(trialId)/(sf*60),pumpBit)
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


