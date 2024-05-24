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
exptName = [animal '_u' unit '_' expt];
exptDir = fullfile(physDir,animal,exptName);
anaMode = 'MU';

load(fullfile(exptDir,[exptName '_id.mat']))
load(fullfile(exptDir,[exptName '_trialInfo.mat']))
if strcmp(anaMode,'SU')
    load(fullfile(exptDir,[exptName '_p' num2str(probe) '_spkSort.mat']))
    spkStruct = spkSort;
    clear spkSort
elseif strcmp(anaMode,'MU')
    load(fullfile(exptDir,[exptName '_p' num2str(probe) '_MUspkMerge.mat']))
    spkStruct = MUspkMerge;
    clear MUspkMerge
end
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
if strcmp(anaMode,'SU')
    nU = length(unique(spkStruct.unitid));
elseif strcmp(anaMode,'MU')
    nU = length(unique(spkStruct.detChSort));
end
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

for u = 4
    if strcmp(anaMode,'SU')
        uIdx = spkStruct.unitid==u;
    elseif strcmp(anaMode,'MU')
        uIdx = spkStruct.detChSort==u;
    end

    win = 0.01;%ms
    nWin = trialL/win;
    for t = 1:nTrials
    
        [rep,cnd] = find(sortTrialInd==t);
        tv = trialStart(t):trialStart(t)+(trialL*sf)-1;
        s(t,:) = zeros(size(tv));
        if cnd ~= trialInfo.blankId
            s(t,ismember(tv,stimStart(t):stimEnd(t))) =  cnd;
        end
        spkTimes{t} = spkStruct.spktimes( uIdx & ( ...
                                        spkStruct.spktimes>trialStart(t) & ...
                                        spkStruct.spktimes<trialEnd(t)) );
        r(t,:) = ismember(tv,spkTimes{t});
        for w = 1:nWin
            S(t,w) = max(unique(s(t,1+((w-1)*(sf*win)):sf*win)));
    %         R(t,w) = 
        end
    
    
    
    end
    
    figure; 
    for c = 1:nConds
        subplot(nConds,1,c);hold on
        tIdx = trialInfo.triallist == c;
        plot(s(tIdx,:)','LineWidth',2)
        plot(r(tIdx,:)','LineWidth',2)
    
    
    end


end