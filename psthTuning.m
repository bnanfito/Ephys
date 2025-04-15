% psth_tuning

clear all
close all

load('/Volumes/NielsenHome2/Brandon/data/dataSets/training/Train_BiDir/MU/threshold4/ranksum/Train_BiDir_MUdataSet.mat')
dat = data.pssbf;
goodIdx = screenUnits(dat,'MU','ranksum');
dat = dat(goodIdx,:);

nU = height(dat);
binSize = 0.1;
bins = -1:binSize:2;
for u = 1:nU
    disp(u)
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
    shift = nConds-find(abs(cVec)==180);
    cVecShift = circshift(cVec,shift);
    cVecShift(cVecShift>180) = cVecShift(cVecShift>180)-360;
    cVecShift(cVecShift<-180) = cVecShift(cVecShift<-180)+360;
    cVecShift = [-180;cVecShift];
    psthShift = circshift(psth,shift);
    psthShift = [psthShift(end,:);psthShift];

%     figure;
%     imagesc(psth)
%     yticks([1:nConds])
%     yticklabels(cId)
% 
%     figure;
%     imagesc(psthShift)
%     yticks([1:nConds])
%     yticklabels(cVecShift)

    psthTune(:,:,u) = psthShift/max(psthShift(:));

end

% [x,y] = meshgrid(bins(2:end),cVecShift);
% surf(x,y,mean(psthTune,3))
imagesc(mean(psthTune,3))
