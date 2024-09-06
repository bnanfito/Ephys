%Kalatsky_BN
clear all
close all

% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
% dataFold = '/Volumes/Lab drive/Brandon/data';
dataFold = 'D:\data';
physDir = fullfile(dataFold,'Ephys');

animal = 'febn1';
unit = '001';
expt = '003';
probe = 2;
exptName = [animal '_u' unit '_' expt];
anaMode = 'SU';

load(fullfile(physDir,animal,exptName,[exptName '_trialInfo.mat']))
load(fullfile(physDir,animal,exptName,[exptName '_id.mat']))
if strcmp(anaMode,'SU')
    load(fullfile(physDir,animal,exptName,[exptName '_p' num2str(probe) '_spkSort.mat']))
elseif strcmp(anaMode,'MU')
    load(fullfile(physDir,animal,exptName,[exptName '_p' num2str(probe) '_MUspkMerge.mat']))
end
load(fullfile(physDir,animal,exptName,[exptName '.analyzer']),'-mat')

%% compute event times

predelay = getparam('predelay',Analyzer);
stimTime = getparam('stim_time',Analyzer);
postdelay = getparam('postdelay',Analyzer);
trialL = sum([predelay,stimTime,postdelay]);

nTrials = length(trialInfo.triallist);
nConds = length(unique(trialInfo.triallist));
[trialSort,trialSortInd] = sort(trialInfo.triallist);
tv = (0:stimTime*id.sampleFreq)/id.sampleFreq;
if strcmp(anaMode,'SU')
    nUnits = length(spkSort.unitinfo);
elseif strcmp(anaMode,'MU')
    nUnits = id.probes(probe).nChannels;
end
kernel = normpdf(-3:6/2000:3);

eventTimes = trialInfo.eventTimes;
eventIDs = trialInfo.eventId;
d = diff(eventIDs);
trialStartInd = [1;find(d==1)+1];
stimStartInd = [2;find(d==1)+2];
stimEndInd = find(d==-1)-1;
trialEndInd = find(d==-1)+1;

eventInds = [trialStartInd stimStartInd stimEndInd trialEndInd];

for t = 1:nTrials

    trialStimInds = ismember(1:length(eventIDs),eventInds(t,2):eventInds(t,3))';

    cycleStart{t} = eventTimes(eventIDs==3 & trialStimInds);
    cycleEnd{t} = cycleStart{t}+round(mean(diff(cycleStart{t})));
    
end

for c = 1:nConds

    cyclePeriod(c) = round(mean( vertcat(cycleEnd{trialInfo.triallist==c}) - vertcat(cycleStart{trialInfo.triallist==c}) ));
    cycleFreq(c) = 1/(cyclePeriod(c)/id.sampleFreq);
    harmonics(:,c) = cycleFreq(c)*(1);

end


%% compute sdf

for u = 1:nUnits

    for t = 1:nTrials
    
        tStart = eventTimes(eventInds(t,1));
        sStart = eventTimes(eventInds(t,2));
        sEnd = eventTimes(eventInds(t,3));
        tEnd = eventTimes(eventInds(t,4));
    
        trialTV = tStart:tStart+(trialL*id.sampleFreq);
        stimTV = sStart:sStart+(stimTime*id.sampleFreq);
    
        if strcmp(anaMode,'SU')
            spkTrain(:,t,u) = ismember(stimTV,spkSort.spktimes(spkSort.unitid==u));
        elseif strcmp(anaMode,'MU')
            spkTrain(:,t,u) = ismember(stimTV,MUspkMerge.spktimes(MUspkMerge.detCh==u));
        end

        for cycle = 1:length(cycleStart{t})

            if isnan(cycleEnd{t}(cycle))
                continue
            end

            cycleTV = cycleStart{t}(cycle):cycleStart{t}(cycle)+cyclePeriod(trialInfo.triallist(t));
            if strcmp(anaMode,'SU')
                spkTrainCycle{t,u}(:,cycle) = ismember(cycleTV,spkSort.spktimes(spkSort.unitid==u));
            elseif strcmp(anaMode,'MU')
                spkTrainCycle{t,u}(:,cycle) = ismember(cycleTV,MUspkMerge.spktimes(MUspkMerge.detCh==u));
            end
            if trialInfo.triallist(t) == 3 || trialInfo.triallist(t) == 4
                spkPhase{t,cycle,u} = ( find(  spkTrainCycle{t,u}(:,cycle)  )/cyclePeriod(trialInfo.triallist(t)) )*360;
            elseif trialInfo.triallist(t) == 1 || trialInfo.triallist(t) == 2
                spkPhase{t,cycle,u} = ( find(  flipud(spkTrainCycle{t,u}(:,cycle))  )/cyclePeriod(trialInfo.triallist(t)) )*360;
            end

        end
    
    end
    
    for c = 1:nConds

        spkTrainCond(:,c,u) = sum(spkTrain(:,trialInfo.triallist==c,u),2);
        sdfCond(:,c,u) = conv(spkTrainCond(:,c,u),kernel,'same');
    
    end

end

for u = 1:nUnits
    figure;
    subplot(3,2,1);hold on
    [x,y] = find(spkTrain(:,trialSortInd,u));
    plot(x,y,'.')

    for c = 1:nConds
        if c==trialInfo.blankId
            continue
        end
        condSpkTrainMat = [spkTrainCycle{trialInfo.triallist==c,u}];

        if c == 1 || c == 3
            col = 'r';
        elseif c == 2 || c == 4
            col = 'b';
        end

        subplot(12,2,2*c);hold on
        [x,y] = find(condSpkTrainMat);
        plot(x,y,[col '.'])

        sdfCycle{c}(:,u) = conv(sum(condSpkTrainMat,2),kernel,'same');
        plot((sdfCycle{c}(:,u)/max(sdfCycle{c}(:,u)))*max(y),col,'LineWidth',2)
    end

%     subplot(3,2,3);hold on
%     for c = 1:nConds
%         plot(tv,sdfCond(:,c,u),'LineWidth',2)
%     end
%     legend({'cond 1','cond 2','cond 3','cond 4'})
    
    subplot(3,2,4);hold on

    yPer = min(cyclePeriod([2,4]));
    yPos = 1:yPer;
    xPer = min(cyclePeriod([1,3]));
    xPos = 1:xPer;

    sdfUpDown = sdfCycle{4}(yPos,u) + flipud(sdfCycle{2}(yPos,u));
    sdfLeftRight = sdfCycle{3}(xPos,u) + flipud(sdfCycle{1}(xPos,u));
    updownPhaseDist = vertcat(spkPhase{trialInfo.triallist==2|trialInfo.triallist==4,:,u});
    leftrightPhaseDist = vertcat(spkPhase{trialInfo.triallist==1|trialInfo.triallist==3,:,u});
    gY = gauss(0:360,mean(updownPhaseDist),std(updownPhaseDist));
    gX = gauss(0:360,mean(leftrightPhaseDist),std(leftrightPhaseDist));

    histogram( updownPhaseDist, 360/2 ,'FaceColor','b','FaceAlpha',0.5,'EdgeColor','none')
    plot(0:360,gY*max(sdfUpDown),'--','LineWidth',2,'Color','b');
    plot((yPos/max(yPos))*360,sdfUpDown,'LineWidth',2,'Color','b')

    histogram( leftrightPhaseDist, 360/2 ,'FaceColor','r','FaceAlpha',0.5,'EdgeColor','none')
    plot(0:360,gX*max(sdfLeftRight),'--','LineWidth',2,'Color','r');
    plot((xPos/max(xPos))*360,sdfLeftRight,'LineWidth',2,'Color','r')

    xlim([0 360])

    subplot(3,2,3);hold on
    for c = 1:nConds
        [X,f,P,Nt,phase]=myfft(sdfCond(:,c,u),id.sampleFreq,0);
        plot(f(f<10),P(f<10))
        for h = 1:size(harmonics,1)
            H(h) = f(abs(f-harmonics(h,c)) == min( abs(f-harmonics(h,c)) ));
        end
        plot(f(ismember(f,H)),P(ismember(f,H)),'o')
        Pow(c,:) = P(ismember(f,H));
        Pha(c,:) = mod(phase(ismember(f,H)),360);
        if c == 1 || c == 2
           Pha(c,:) = abs(Pha(c,:)-360)+1;
        end
    end

    subplot(3,2,5);hold on
    x_pha = Pha([1 3],:); x_pha = x_pha(:);
    y_pha = Pha([2 4],:); y_pha = y_pha(:);
    x_pow = Pow([1 3],:); x_pow = x_pow(:);
    y_pow = Pow([2 4],:); y_pow = y_pow(:);
    wmeanX_pha = sum(x_pha.*x_pow)/sum(x_pow);
    wmeanY_pha = sum(y_pha.*y_pow)/sum(y_pow);
    
    gX_pha = gauss(0:360,wmeanX_pha,std(x_pha));
    gY_pha = gauss(0:360,wmeanY_pha,std(y_pha));
    imagesc(gY_pha'*gX_pha)
    xlim([0 360])
    ylim([0 360])

    subplot(3,2,6)
%     imagesc(gY'*gX)
%     xlim([0 360])
%     ylim([0 360])
    
    meshSDF(:,:,u) =  downsample(sdfUpDown,10)*downsample(sdfLeftRight,10)';
    meshSDF_x(:,u) = (downsample(xPos,10)/max(downsample(xPos,10)))*360;
    meshSDF_y(:,u) = (downsample(yPos,10)/max(downsample(yPos,10)))*360;
    imagesc(meshSDF_x(:,u),meshSDF_y(:,u),meshSDF(:,:,u))




end

figure;imagesc(mean(meshSDF_x,2),mean(meshSDF_y,2),mean(meshSDF,3))
