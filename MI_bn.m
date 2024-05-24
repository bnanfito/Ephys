% Mutual information analysis
clear all
close all

if ispc
%     dataFold = 'D:\data'; 
%     dataFold = 'C:\Users\brand\Documents\data';
    dataFold = 'F:\Brandon\data';
elseif ismac
%     dataFold = '/Volumes/Lab drive/Brandon/data';
%     dataFold = '/Volumes/Lab drive/Brandon/VSS2024/data';
    dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end
physDir = fullfile(dataFold,'Ephys');

animal = 'febl0';
unit = '000';
expt = '010';
probe = 1;
u = 1;
exptName = [animal '_u' unit '_' expt];
exptDir = fullfile(physDir,animal,exptName);

load(fullfile(exptDir,[exptName '_id.mat']))
load(fullfile(exptDir,[exptName '_trialInfo.mat']))
load(fullfile(exptDir,[exptName '_p' num2str(probe) '_spkSort.mat']))
load(fullfile(exptDir,[exptName '.analyzer']),'-mat')

predelay = getparam('predelay',Analyzer);
stimTime = getparam('stim_time',Analyzer);
postdelay = getparam('postdelay',Analyzer);
trialL = predelay+stimTime+postdelay;

sf = id.sampleFreq;
nTrials = length(trialInfo.triallist);
nConds = length(unique(trialInfo.triallist));
nReps = nTrials/nConds;
nEpochs = length(trialInfo.eventTimes)/nTrials;
if isempty(trialInfo.blankId)
    blank = zeros(1,nConds)==1;
else
    blank = (1:nConds)==trialInfo.blankId;
end

eventID = [1;diff(trialInfo.eventId)];
eventTimes = trialInfo.eventTimes;
trialStart = eventTimes(eventID == 1);
stimStart = eventTimes(eventID == 2);
stimEnd = eventTimes(eventID == -2);
trialEnd = eventTimes(eventID == -1);

[sortTrialCond,sortTrialInd] = sort(trialInfo.triallist);
sortTrialCond = reshape(sortTrialCond,nReps,nConds);
sortTrialInd = reshape(sortTrialInd,nReps,nConds);

win = 0.01;%ms
nWin = trialL/win;
for t = 1:nTrials

    [r,c] = find(sortTrialInd==t);
    trialInfo.triallist(t) == c
    tv = trialStart(t):trialStart(t)+(trialL*sf);
    s(t,:) = zeros(size(tv));
    if trialInfo.triallist(t) ~= trialInfo.blankId
        s(ismember(tv,stimStart(t):stimEnd(t))) =  trialInfo.triallist(t);
    end
    spkTimes{t} = spkSort.spktimes( spkSort.unitid==u & ( ...
                                    spkSort.spktimes>trialStart(t) & ...
                                    spkSort.spktimes<trialEnd(t)) );
    r(t,:) = ismember(tv,spkTimes{t});
    for w = 1:nWin
        S(t,w) = max(unique(s(t,1+((w-1)*(sf*win)):sf*win)));
%         R(t,w) = 
    end



end



