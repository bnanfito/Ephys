% psth_tuning

clear all
close all

load('Y:\Brandon\data\dataSets\training\Train_V1Cool\MU\threshold4\ranksum\Train_V1Cool_MUdataSet.mat')
dat = data.pssaf;
goodIdx = screenUnits(dat,'MU','ranksum');
dat = dat(goodIdx,:);

nU = height(dat);
binSize = 0.1;
bins = -1:binSize:2;
for u = 1:nU
    [cndInclude] = limitConds(dat(u,:));
    nConds = sum(cndInclude);
    nReps = size(dat.fr(u).trialNum,1);
    nTrials = nReps*nConds;
    conds = unique(dat.fr(u).trialCond(:,cndInclude));
    trials = dat.fr(u).trialNum(:,cndInclude); trials = trials(:);
    spkTrials = dat.spkTimes{u}(2,:);
    spkTimes = dat.spkTimes{u}(1,:);
    
    for c = 1:nConds
        % list of trial numbers that are of the current condition
        condTrials = dat.fr(u).trialNum(dat.fr(u).trialCond == conds(c));
        condTrials = condTrials(:);
        condIdx = ismember(spkTrials,condTrials);

        psth(c,:) = histcounts(spkTimes(condIdx),bins);
        psth(c,:) = psth(c,:)/(nReps*binSize);
    end
    cId = dat.cndKey{u}(cndInclude,strcmp(dat.paramKey{u},'ori'));
    cPref = dat.oriPref(u);
    cVec = cId-cPref;
    cVec(cVec<-180) = cVec(cVec<-180)+360;
    cVec(cVec>180) = cVec(cVec>180)-360;
    [cVecShift,shiftIdx] = sort(cVec);
    psthShift = psth(shiftIdx,:);
    if abs(cVecShift(1)) < abs(cVecShift(end))
        cVecShift = [cVecShift(end)*-1;cVecShift];
        shiftIdx = [shiftIdx(end);shiftIdx];
        psthShift = [psthShift(end,:);psthShift];
    elseif abs(cVecShift(1)) > abs(cVecShift(end))
        cVecShift = [cVecShift;cVecShift(1)*-1];
        shiftIdx = [shiftIdx;shiftIdx(1)];
        psthShift = [psthShift;psthShift(1,:)];
    end

%     figure;
%     imagesc(psth)
%     yticks([1:nConds])
%     yticklabels(cId)
% 
%     figure;
%     imagesc(psthShift)
%     yticks([1:nConds+1])
%     yticklabels(cVecShift)

    psthTune(:,:,u) = psthShift/max(psthShift(:));

end

% [x,y] = meshgrid(bins(2:end),cVecShift);
% surf(x,y,mean(psthTune,3))
imagesc(mean(psthTune,3))

psthPref = squeeze(psthTune(cVecShift == 0,:,:))';
psthNull = squeeze(psthTune(cVecShift == 180,:,:))';
psthOrth1 = squeeze(psthTune(cVecShift == 90,:,:))';
psthOrth2 = squeeze(psthTune(cVecShift == -90,:,:))';

figure; hold on
plot(mean(psthPref))
plot(mean(psthNull))
plot(mean(psthOrth1))
plot(mean(psthOrth2))

figure; hold on
plot(mean(psthTune,3)')