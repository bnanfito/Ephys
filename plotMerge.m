%plot merged
clear all
close all

dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
animal = 'febl0';
units = {'000','000'};
expts = {'010','012'};
probe = 1;

physDir = fullfile(dataFold,'Ephys');
clr = {'k','b'};
for f = 1:length(expts)

    exptName = [animal '_u' units{f} '_' expts{f}];
    load(fullfile(physDir,animal,exptName,[exptName '_id.mat']),'id')
    load(fullfile(physDir,animal,exptName,[exptName '_trialInfo.mat']),'trialInfo')
    load(fullfile(physDir,animal,exptName,[exptName '.analyzer']),'-mat')
    load(fullfile(physDir,animal,exptName,[animal '_u' units{f} '_' expts{f} '_p' num2str(probe) '_' num2str(f) '_spkSort.mat']),'spkSort')
    
    area = id.probes(probe).area;
    sf = id.sampleFreq;
    % kernel = ones(1,0.1*sf)*(1/(0.1*sf));
    kernel = normpdf(-3:6/2000:3);
    colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
    nTrials = length(trialInfo.triallist);
    nConds = length(unique(trialInfo.triallist));
    nDom = length(trialInfo.dom);
    nReps = nTrials/nConds;
    nU = length(spkSort.unitinfo);
    uIDs = unique(spkSort.unitid); uIDs = uIDs(uIDs~=0);

    predelay = getparam('predelay',Analyzer);
    stimTime = getparam('stim_time',Analyzer);
    postdelay = getparam('postdelay',Analyzer);
    if stimTime < 1
        st = stimTime;
    else
        st = 1;
    end
    trialStart = trialInfo.eventTimes([1;find(diff(trialInfo.eventId) == 1)+1]);
    stimStart = trialInfo.eventTimes(find(diff(trialInfo.eventId) == 2)+1);
    stimEnd = trialInfo.eventTimes(find(diff(trialInfo.eventId) == -2)+1);
    trialEnd = trialInfo.eventTimes(find(diff(trialInfo.eventId) == -1)+1);
    [sortTrialCond,sortTrialInd] = sort(trialInfo.triallist);
    sortTrialCond = reshape(sortTrialCond,nReps,nConds);
    sortTrialInd = reshape(sortTrialInd,nReps,nConds);
    trialL = predelay+stimTime+postdelay;
    nT = ((1:(trialL*sf))-(predelay*sf))/sf;
    if isempty(trialInfo.blankId)
        blank = zeros(1,nConds)==1;
    else
        blank = (1:nConds)==trialInfo.blankId;
    end

    for u = uIDs
        figure(u);hold on
        spkTimes = spkSort.spktimes(spkSort.unitid == u);
        for c = 1:nConds
            trials = find(trialInfo.triallist == c);
            for r = 1:nReps
                t = trials(r);
                tvPre       = stimStart(t)-(predelay*sf):stimStart(t)-1;
                tvStim      = stimStart(t):stimStart(t)+(stimTime*sf)-1;
                tvPost      = stimStart(t)+(stimTime*sf):stimStart(t)+((stimTime+postdelay)*sf)-1;
                tvTrial     = [tvPre tvStim tvPost];
                spkTrain{f}(:,u) = ismember(tvTrial,spkTimes);
                baseCount = length(find(ismember(stimStart(t)-sf:stimStart(t),spkTimes)));
                stimCount = length(find(ismember(stimStart(t):stimStart(t)+sf,spkTimes)));
                baseFR{f}(r,c,u) = baseCount;
                stimFR{f}(r,c,u) = stimCount;
                bcfr{f}(r,c,u) = stimCount-baseCount;

            end

        end
        x = trialInfo.domval(sortTrialCond(:,~blank));
        y = bcfr{f}(:,:,u);
        plot(x,y(:,~blank),[clr{f} '.'])
        patch([x(1,1) x(1,end) x(1,end) x(1,1)],[min(y(:,blank)) min(y(:,end)) max(y(:,end)) max(y(:,end))],'r','EdgeColor','none','FaceAlpha',0.2)
        plot(mean(x),mean(y(:,~blank)),[clr{f} '-o'])



    end


end