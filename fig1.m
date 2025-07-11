clear all
close all

animal = 'febl4';
unit = '000';
expt = '012';
probe = 1;
anaMode = 'MU';
dataFold = 'F:\Brandon\data';
exptName = [animal '_u' unit '_' expt];

% load necessary data
load(fullfile(dataFold,'Ephys',animal,exptName,[exptName '_temp.mat']))
load(fullfile(dataFold,'Ephys',animal,exptName,[exptName '_trialInfo.mat']))
load(fullfile(dataFold,'Ephys',animal,exptName,[exptName '_id.mat']))
load(fullfile(dataFold,'Ephys',animal,exptName,[exptName '.analyzer']),'-mat')
load('F:\Brandon\data\Ephys\febl4\febl4_u000_011\goodUnits.mat')
[sumStats,spks] = anaOri('febl4','000','012',1,'MU',dataFold,0,0);

pre = getparam('predelay',Analyzer);
stimTime = getparam('stim_time',Analyzer);
post = getparam('postdelay',Analyzer);

sf = id.sampleFreq;
nTrials = length(trialInfo.triallist);
goodIdx = ismember([spks.unitId],goodUnitId);
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

% for u = 1:height(sumStats)
%     coldR(:,u) = sumStats.fr(u).bc(ismember(sumStats.fr(u).trialNum,coldTrials));
%     warmR(:,u) = sumStats.fr(u).bc(ismember(sumStats.fr(u).trialNum,warmTrials));
% end
% warmR = mean(warmR(:,goodIdx),1,'omitnan');
% warmR(warmR<0) = 0;
% maxR = max(warmR);
% warmR = warmR/maxR;
% warmR = 1-warmR;
% coldR = mean(coldR(:,goodIdx),1,'omitnan');
% coldR(coldR<0) = 0;
% coldR = coldR/maxR;
% coldR = 1-coldR;

eventIDdiff = [1;diff(trialInfo.eventId)];
eventTimes = trialInfo.eventTimes;
trialStart = eventTimes(eventIDdiff==1);
stimStart = eventTimes(eventIDdiff==2);
stimEnd = eventTimes(eventIDdiff==-2);
trialEnd = eventTimes(eventIDdiff==-1);

%% plot

for u = goodUnitId'

figure;
tiledlayout(2,2)
bw = id.sampleFreq*1;

uIdx = [spks.unitId] == u;
nexttile; hold on
spkTimes = spks(uIdx).stimCent;
patch([0 1 1 0],[0 0 nTrials+1 nTrials+1],'k','EdgeColor','none','FaceAlpha',0.2)
scatter(spkTimes(1,:),spkTimes(2,:),1,'k.')
ylim([0 nTrials+1])
ylabel('trial')
xlim([-1 2])
xlabel('time (sec)')
title(['channel ' num2str(spks(uIdx).unitId)])
box on

nexttile; hold on
colororder({'k','r'})
[psth,edges] = histcounts(spks(uIdx).times,'BinWidth',bw);
yyaxis left
plot(edges(2:end)/(sf*60),psth)
ylabel('firing rate (hz)')
yyaxis right
plot(trialStart(trialId)/(sf*60),tempDat)
ylabel('temperature (c)')
xlabel('time (min)')
box on

nexttile; hold on
spkTimes = [spks(:).stimCent];
patch([0 1 1 0],[0 0 nTrials+1 nTrials+1],'k','EdgeColor','none','FaceAlpha',0.2)
scatter(spkTimes(1,:),spkTimes(2,:),1,'k.')
ylim([0 nTrials+1])
ylabel('trial')
xlim([-1 2])
xlabel('time (sec)')
title("'good' channels")
box on

nexttile; hold on
[psth,edges] = histcounts([spks(goodIdx).times],'BinWidth',bw);
psth = psth/sum(goodIdx);
yyaxis left
plot(edges(2:end)/(sf*60),psth)
ylabel('firing rate (hz)')
yyaxis right
plot(trialStart(trialId)/(sf*60),tempDat)
ylabel('')
xlabel('time (min)')
box on

sgtitle([animal ' ' unit ' ' expt])

end

% figure; hold on
% scatter([spks.xPos],[spks.zPos])
% scatter([spks(~goodIdx).xPos],[spks(~goodIdx).zPos],'rx')
% scatter(sumStats.xPos(goodIdx),sumStats.zPos(goodIdx),200,repmat( warmR', 1,3 ),'.')
% 
% figure; hold on
% scatter([spks.xPos],[spks.zPos])
% scatter([spks(~goodIdx).xPos],[spks(~goodIdx).zPos],'rx')
% scatter(sumStats.xPos(goodIdx),sumStats.zPos(goodIdx),200,repmat( coldR', 1,3 ),'.')


