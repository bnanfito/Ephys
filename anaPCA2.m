% clear all
% close all

function [r,cid,tid,bid,binEdges,acc,dat] = anaPCA2(dat)


bw = .25;
tBase = 1;
tStim = 2;
binEdges = -tBase:bw:tStim;
nBins = length(binEdges)-1;

% load('V1cool_SU_ori_projectTbl.mat')
% coolIdx = projTbl.duringMFlag == 1;
% ageIdx = projTbl.age>=34 & projTbl.age<=37;
% % dat = vertcat(projTbl.sumStats{61});
% dat = vertcat(projTbl.sumStats{~coolIdx & ageIdx});

load('DSdev_dataSet.mat')
cntrlIdx = projectTbl.duringMFlag == 0 & projectTbl.priorMFlag == 0;
ageIdx = projectTbl.age>=34 & projectTbl.age<=37;
areaIdx = strcmp(projectTbl.recSite,'PSS');
dat = vertcat(projectTbl.sumStats{cntrlIdx & ageIdx & areaIdx});



dat = dat(dat.goodUnit,:);
[~,oriPrefIdx] = sort(dat.oriPref);
dat = dat(oriPrefIdx,:);
nU = height(dat);
for u = 1:nU

    spkTimes = dat.spkTimes{u}(1,:); %vector of spike times for unit u
    spkTrialID = dat.spkTimes{u}(2,:); %trial ID for each entry in spkTimes
    %if stimulus size is a changing parameter, only include stimulus conditions with small size (exclude fullfield conditions)
    if sum(contains(dat.paramKey{u},'size'))==1 
        cndInclude = find(dat.cndKey{u}(:, contains(dat.paramKey{u},'size') )<100)'; 
    else
        cndInclude = 1:length(unique(dat.cndKey{u}(:, contains(dat.paramKey{u},'ori') )));
    end
    cndVal = dat.cndKey{u}(cndInclude, contains(dat.paramKey{u},'ori') );
    trialID = dat.fr(u).trialNum(:,cndInclude);
    condID = repmat(1:length(cndInclude),size(trialID,1),1);
    trialID = trialID(:);
    condID = condID(:);
    nTrials = length(trialID);

    spkInclude = ismember(spkTrialID,trialID); %only include spikes from trials with small stim condition (exclude fullfield trials)
    spkTimes = spkTimes(spkInclude);
    spkTrialID = spkTrialID(spkInclude);

    for t = 1:nTrials
        curSpks = spkTimes(spkTrialID==trialID(t));
        for b = 1:nBins
            rBoxCar(b,t,u) = sum( curSpks>binEdges(b)-(bw) & curSpks<binEdges(b)+(bw) );
        end
        [rHist(:,t,u),~] = histcounts(curSpks,binEdges);
    end
    
    spks{u} = spkTimes;

end

% r = rBoxCar;
r = rHist;
cid = repmat(condID',nBins,1);
tid = repmat(trialID',nBins,1);
bid = repmat((1:nBins)',1,nTrials);

R = reshape(r,nBins*nTrials,nU);
R = R./max(R);
Cid = reshape(cid,nBins*nTrials,1);
Tid = reshape(tid,nBins*nTrials,1);
Bid = reshape(bid,nBins*nTrials,1);
trials = unique(Tid);
[coeff,score,latent,tsquare,explained] = pca(R);

for b = 1:size(r,1)

    [acc(b),~] = popDecode(squeeze(r(b,:,:)),cid(b,:));

end



%% plot

figure;hold on
histogram([spks{:}],binEdges)
title('PSTH (all trials/units)')
ylabel('spike count')
xlabel('time (s) rel. to stim onset')

figure;hold on
pcInclude = [2 3 4];
conds = unique(Cid);
clrs = hsv(length(conds));
for cond = 1:length(conds)
    p(cond) = plot3(score(Cid==cond,pcInclude(1)),score(Cid==cond,pcInclude(2)),score(Cid==cond,pcInclude(3)),'.','Color',clrs(cond,:),'MarkerSize',10);
end
for t = 1:length(trials)
%     plot3(score(Tid==trials(t),pcInclude(1)),score(Tid==trials(t),pcInclude(2)),score(Tid==trials(t),pcInclude(3)),'k--')
    plot3(score(Tid==trials(t),pcInclude(1)),score(Tid==trials(t),pcInclude(2)),score(Tid==trials(t),pcInclude(3)),'--','Color',clrs(unique(Cid(Tid==trials(t))),:),'LineWidth',1)
end
legend(num2str(cndVal))
xlabel(['PC' num2str(pcInclude(1))])
ylabel(['PC' num2str(pcInclude(2))])
zlabel(['PC' num2str(pcInclude(3))])

figure;hold on
stem(cumsum(explained))
xlabel('PC')
ylabel('cum.sum explained variance')
ylim([0 100])

figure;hold on
imagesc(R)
axis tight
ylabel('trial/time bin')
xlabel('unit')
colorbar

figure;hold on
plot(acc)
ylabel('accuracy')
xlabel('time bin')
title('decoder performance')




end
