% Written by Brandon Nanfito
function [spks] = orgSpks(animal,unit,expt,probe,anaMode,dataFold)

    baseName = [animal '_u' unit '_' expt];
%     physDir = fullfile(dataFold,'Ephys',animal,baseName);
%     physDir = fullfile(dataFold,'sf4rs01',animal,baseName);
    physDir = fullfile(dataFold,'sf5rs01',animal,baseName);
    load(fullfile(physDir,[baseName '_id.mat']),'id');
    load(fullfile(physDir,[baseName '_trialInfo.mat']),'trialInfo');
    load(fullfile(physDir,[baseName '.analyzer']),'-mat','Analyzer');
    if strcmp(anaMode,'SU')
        load(fullfile(physDir,[baseName '_p' num2str(probe) '_spkSort.mat']),'spkSort');
        spkStruct = spkSort; 
        clear spkSort
    elseif strcmp(anaMode,'MU')
        load(fullfile(physDir,[baseName '_p' num2str(probe) '_MUspkMerge.mat']),'MUspkMerge');
        spkStruct = MUspkMerge;
        clear MUspkMerge
    end

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
    partL = 1; %second of stimulus to calculate FR from if stimTimes>2sec (avoid habituated period)

    predelay = getparam('predelay',Analyzer);
    stimTime = getparam('stim_time',Analyzer);
    postdelay = getparam('postdelay',Analyzer);

    trialStart = downsample(trialInfo.eventTimes,nEpochs);
    stimStart = downsample(trialInfo.eventTimes,nEpochs,1);
    stimEnd = downsample(trialInfo.eventTimes,nEpochs,2);
    trialEnd = downsample(trialInfo.eventTimes,nEpochs,3);

    if strcmp(anaMode,'SU')
        nUnits = max(spkStruct.unitid);
    elseif strcmp(anaMode,'MU')
        nUnits = id.probes(probe).nChannels;
    end

    for u = 1:nUnits 
    
        if strcmp(anaMode,'SU')
            spks(u).info = spkStruct.unitinfo{u};
            spks(u).times = spkStruct.spktimes(spkStruct.unitid == u);
        elseif strcmp(anaMode,'MU')
            spks(u).info = 'MU';
            spks(u).times = spkStruct.spktimes(spkStruct.detCh==u);
        end

        spks(u).stimCent = [];

        for c = 1:nConds
            trials = find(trialInfo.triallist == c);
            for r = 1:nReps
                t = trials(r);
                
                % time vectors for different trial epochs
                tvPre       = stimStart(t)-(predelay*sf):stimStart(t)-1;
                tvStim      = stimStart(t):stimStart(t)+(stimTime*sf)-1;
                tvStim_part = stimStart(t)+((partL-1)*sf):stimStart(t)+(partL*sf)-1;
                tvPost      = stimStart(t)+(stimTime*sf):stimStart(t)+((stimTime+postdelay)*sf)-1;
                tvTrial     = [tvPre tvStim tvPost];

                spks(u).train(:,t) = ismember(tvTrial,spks(u).times);
                spks(u).stimCent = [spks(u).stimCent vertcat((tvTrial(spks(u).train(:,t))-stimStart(t))/sf,repmat(t,1,length(find(spks(u).train(:,t)))))];

                baseSpkCount = length(find(ismember(tvPre,spks(u).times)));
                baseFR = baseSpkCount/predelay;
                spks(u).fr.base(r,c) = baseFR;

                if stimTime>2
                    stimSpkCount = sum(ismember(tvStim_part,spks(u).times));
                    stimFR = stimSpkCount/1;
                    spks(u).fr.stim(r,c) = stimFR;
                else
                    stimSpkCount = sum(ismember(tvStim,spks(u).times));
                    stimFR = stimSpkCount/stimTime;
                    spks(u).fr.stim(r,c) = stimFR;
                end

                spks(u).fr.bc(r,c) = stimFR-baseFR;

            end
        end
    
    end


end