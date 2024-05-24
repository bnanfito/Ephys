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
expt = '005';
probe = 1;
exptName = [animal '_u' unit '_' expt];
exptDir = fullfile(physDir,animal,exptName);

load(fullfile(exptDir,[exptName '_id.mat']))
load(fullfile(exptDir,[exptName '_trialInfo.mat']))
load(fullfile(exptDir,[exptName '_p' num2str(probe) '_spkSort.mat']))
load(fullfile(exptDir,[exptName '.analyzer']),'-mat')

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
stimPartL = 1; %if stimTime>1second --> only use stimPartL seconds of the stim
stimPart = 1; %which stimPartL seconds to use; 1 = 1st stimPartL seconds

predelay = getparam('predelay',Analyzer);
stimTime = getparam('stim_time',Analyzer);
postdelay = getparam('postdelay',Analyzer);
if stimTime < 1
    st = stimTime;
else
    st = 1;
end
eventID = sum(trialInfo.eventCh,2);
trialStart = downsample(trialInfo.eventTimes,nEpochs);
stimStart = downsample(trialInfo.eventTimes,nEpochs,1);
stimEnd = downsample(trialInfo.eventTimes,nEpochs,2);
trialEnd = downsample(trialInfo.eventTimes,nEpochs,3);
[sortTrialCond,sortTrialInd] = sort(trialInfo.triallist);
sortTrialCond = reshape(sortTrialCond,nReps,nConds);
sortTrialInd = reshape(sortTrialInd,nReps,nConds);
