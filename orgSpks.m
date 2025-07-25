
% Written by Brandon Nanfito
function [spks,trialExclude] = orgSpks(animal,unit,expt,probe,anaMode,dataFold)

    baseName = [animal '_u' unit '_' expt];
    physDir = fullfile(dataFold,'Ephys',animal,baseName);
    load(fullfile(physDir,[baseName '_id.mat']),'id');
    load(fullfile(physDir,[baseName '_trialInfo.mat']),'trialInfo');
    load(fullfile(physDir,[baseName '.analyzer']),'-mat','Analyzer');

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

    eventIDdiff = [1;diff(trialInfo.eventId)];
    eventTimes = trialInfo.eventTimes;
    if strcmp(Analyzer.modID,'PG')

        stimStartID = 3;
        trialStart = eventTimes(eventIDdiff==1);
        stimStart = eventTimes(eventIDdiff==2);
        stimEnd = eventTimes(eventIDdiff==-2);
        trialEnd = eventTimes(eventIDdiff==-1);

    elseif strcmp(Analyzer.modID,'RD')

        stimStartID = 7;
        trialStart = eventTimes(eventIDdiff==1);
        stimStart = eventTimes(eventIDdiff==6 | eventIDdiff==2);
        stimEnd = eventTimes(eventIDdiff==-6 | eventIDdiff==-2);
        trialEnd = eventTimes(eventIDdiff==-1);

    elseif strcmp(Analyzer.modID,'BK')
        
%         stimStartID = 1;
%         trialStart = eventTimes(eventIDdiff==1);
%         stimStart = eventTimes(eventIDdiff==2);
%         trialEnd = eventIDdiff(eventIDdiff==-1);

    end
    [sortTrialCond,sortTrialInd] = sort(trialInfo.triallist);
    sortTrialCond = reshape(sortTrialCond,nReps,nConds);
    sortTrialInd = reshape(sortTrialInd,nReps,nConds);


    if strcmp(anaMode,'SU')
        SUTrialData(fullfile(dataFold,'Ephys'),animal,unit,expt,probe,'id',stimStartID,st,st,0)
        load(fullfile(physDir,[baseName '_p' num2str(probe) '_SUTrial.mat']),'SU','SUinfo')
        nUnits = length(SU);
        trialExclude = zeros(nTrials,1) == 1;
        load(fullfile(physDir,[baseName '_p' num2str(probe) '_spkSort.mat']),'spkSort');
        spkStruct = spkSort; 
        latency = computeLatency_bn(SU,SUinfo,trialExclude);
        spkTrialDat = SU;
        clear spkSort SUlatency
    elseif strcmp(anaMode,'MU')
        MUThreshTrialData(fullfile(dataFold,'Ephys'),animal,unit,expt,probe,'id',stimStartID,st,st,0)
        load(fullfile(physDir,[baseName '_p' num2str(probe) '_MUThreshTrial.mat']),'MUThresh','MUThreshInfo')
        nUnits = length(MUThresh);
        trialExclude = MUThreshInfo.trialExclude;
        load(fullfile(physDir,[baseName '_p' num2str(probe) '_MUspkMerge.mat']),'MUspkMerge');
        spkStruct = MUspkMerge;
        latency = computeLatency_bn(MUThresh,MUThreshInfo,trialExclude);
        spkTrialDat = MUThresh;
        clear MUspkMerge MUlatency
    end


    for u = 1:nUnits 
    
        if strcmp(anaMode,'SU')
            spks(u).unitId = spkTrialDat(u).unitId;
            spks(u).info = spkTrialDat(u).unitClass;
            spks(u).times = spkStruct.spktimes(spkStruct.unitid == spks(u).unitId);
        elseif strcmp(anaMode,'MU')
            spks(u).unitId = spkTrialDat(u).detCh;
            spks(u).info = 'MU';
            spks(u).times = spkStruct.spktimes(spkStruct.detCh== spks(u).unitId);
            spks(u).xPos = id.probes(probe).x(spks(u).unitId);
            spks(u).zPos = id.probes(probe).z(spks(u).unitId);
            spks(u).shaft = id.probes(probe).shaft(spks(u).unitId);
        end

        spks(u).stimCent = [];
        spks(u).late = latency(u)/1000;
        spks(u).fr.trialNum = sortTrialInd;
        spks(u).fr.trialCond = sortTrialCond;

        for c = 1:nConds
            trials = find(trialInfo.triallist == c);
            for r = 1:nReps
                t = trials(r);
                
                % time vectors for different trial epochs
                tvPre       = stimStart(t)-(predelay*sf):stimStart(t)-1;
                tvStim      = stimStart(t):stimStart(t)+(stimTime*sf)-1;
                tvPost      = stimStart(t)+(stimTime*sf):stimStart(t)+((stimTime+postdelay)*sf)-1;
                tvTrial     = [tvPre tvStim tvPost];

                spks(u).train(:,t) = ismember(tvTrial,spks(u).times);
                spks(u).stimCent = [spks(u).stimCent vertcat( (tvTrial( spks(u).train(:,t))-stimStart(t))/sf , ...
                                                                repmat(t,1,length(find(spks(u).train(:,t)))) ) ];

                if strcmp(anaMode,'MU') && ismember(t,find(trialExclude))
                    spks(u).fr.base(r,c) = nan;
                    spks(u).fr.stim(r,c) = nan;
                    spks(u).fr.bc(r,c) = nan;
                else
%                     spks(u).fr.base(r,c) = mean(spkTrialDat(u).baseFrate,'omitnan');
                    spks(u).fr.base(r,c) = spkTrialDat(u).baseFrate(t);
                    spks(u).fr.stim(r,c) = spkTrialDat(u).stimFrate(t);
                    spks(u).fr.bc(r,c) = spks(u).fr.stim(r,c)-spks(u).fr.base(r,c);
                end



%                 s = spks(u).stimCent(1,:);
%                 tID = spks(u).stimCent(2,:);
%                 [values,edges] = histcounts(s(s<=1 & tID==t),-predelay:0.05:1);
% %                 values = zscore(values);
%                 values = values(end-19:end);
%                 edges = edges(end-19:end);
%                 spks(u).psth_stim.values(r,:,c) = values;
%                 spks(u).psth_stim.binEdges(r,:,c) = edges;
%                 clear s tID values edges


            end
        end

        meanBase = mean(spks(u).fr.base(:),'omitnan');
        stdBase = std(spks(u).fr.base(:),'omitnan');
        spks(u).fr.z = (spks(u).fr.stim - meanBase)./stdBase;
    
    end

    % remove outlier trials
%     outProbe = sum(cat(3,vertcat(spks(:).fr).out),3)>32;
%     for u = 1:nUnits
%         if strcmp(anaMode,'SU')
%             spks(u).fr.bc(spks(u).fr.out) = nan;
%         elseif strcmp(anaMode,'MU')
%             spks(u).fr.bc(spks(u).fr.out) = nan;
%             spks(u).fr.bc(outProbe) = nan;
%         end
%     end

end