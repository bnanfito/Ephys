%spkCrossCorr.mat - written by Brandon Nanfito

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% animal = animal experiment ID (ex: 'febg0')
% unit = experiment unit number (ex: '000')
% expt = experiment number (ex: '000')
% areaA = string with area name (ex: 'V1') 
% areaB = string with area name (ex: 'PSS')
% plt = binary (1 for yes; 0 for no) indication to plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% peaksAB = matrix storing the peak = max(mean(ccgAB,2)-mean(mean(ccgShAB,3),2)) 
        % in 1st, time lag of peak in 2nd, and difference in stimulus orientation 
        % pref in 3rd position of 3rd dim for each neuron pair; size: number of neurons 
        % in areaA (rows) x number of neurons in areaB (col.s) x 3;
% ccgAB = cell array containing the cross-corellograms for each neuron pair;
        % size: number of neurons in areaA (rows) x number of neurons in areaB (col.s)
        % each cell contains a matrix of size: time lags (rows) x stimulus
        % conditions (col.s)
% ccgShAB = cell array containing the cross-corellograms for each neuron
        % pair using shuffled repeats (within conditions); size: number of
        % neurons in areaA (rows) x number of neurons in areaB (col.s);
        % each cell contains a matrix of size: time lags (rows) x stimulus
        % conditions (col.s) x number of shuffles
% ccgAA = same as ccgAB but with areaA spk trains vs areaA spk trains
% ccgBB = same as ccgAA but with areaB spk trains



function [peaksAB,ccgAB,ccgShAB,ccgAA,ccgBB,rAB,FR,goodUnits,prefCondA,prefCondB] = spkCrossCorr2(animal,unit,expt,areaA,areaB,plt)

anaMode = 'MU';
baseName = [animal,'_u',unit,'_',expt];
if ispc
    dataFold = 'D:\data'; 
%     dataFold = 'C:\Users\brand\Documents\data';
%     dataFold = 'F:\Brandon\data';
elseif ismac
    dataFold = '/Volumes/Lab drive/Brandon/data';
%     dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end
physDir = fullfile(dataFold,'Ephys');
figDir = fullfile(dataFold,'Figures');

% load experiment files
load(fullfile(physDir,animal,baseName,[baseName '_id.mat']))
load(fullfile(physDir,animal,baseName,[baseName '_trialInfo.mat']))
load(fullfile(physDir,animal,baseName,[baseName '.analyzer']),'-mat')

% extract experiment parameters into variables in the workspace
predelay = getparam('predelay',Analyzer);
stimTime = getparam('stim_time',Analyzer);
postdelay = getparam('postdelay',Analyzer);
sf = id.sampleFreq;
nProbes = length(id.probes);
nTrials = length(trialInfo.triallist);
nConds = length(unique(trialInfo.triallist));
nReps = nTrials/nConds;
if isempty(trialInfo.blankId)
    blank = zeros(1,nConds)==1;
else
    blank = (1:nConds)==trialInfo.blankId;
end


% settings
kernel = [0.05 0.25 0.40 0.25 0.05];
smooth = 1;
xCorrMethod = 1;
normMethod = 2; % 1 = overlap * mean FR; 2 = mean FR; 3 = no normalization
binSizeSec = 0.001;
binSizeSam = binSizeSec*sf; % # of samples in each time bin (across which spks will be binned)
timeLengthSec = 0.6; % the length of time from the end of the stim period to be analyzed (in seconds)
timeLengthSam = timeLengthSec*sf;
nBins = floor(timeLengthSam)/binSizeSam; % # of bins (want to analyze 'sf*tl' # of samples in stim period put into '(sf*tl)/bs' # of bins)
lag = -nBins+1:nBins-1;
colors = {[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560]...
         ,[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
colors = repmat(colors,1,10);

% organize event timestamps (samples) and trial info
trialStart = downsample(trialInfo.eventTimes,4);
stimStart = downsample(trialInfo.eventTimes,4,1);
stimEnd = downsample(trialInfo.eventTimes,4,2);
trialEnd = downsample(trialInfo.eventTimes,4,3);
trialConds = trialInfo.triallist;
[sortedTrials,sortTrialInds] = sort(trialConds);
for cond = 1:nConds % identify samples for each trial (indexed by repeat # (row) and stimulus condition (column) in matrix 'nS') 
    condTrials(:,cond) = find(trialConds == cond);
    for rep = 1:nReps

        trial = condTrials(rep,cond);
        nS(:,rep,cond) = stimEnd(trial)-timeLengthSam+1:stimEnd(trial); % time vector (in samples) for creating binary spike vectors

    end
end
nT = ((0:size(nS,1)-1)/sf)+(stimTime-timeLengthSec); % generic time vectors (for plotting purposes)


%% Make spike train matrix

for p = 1:nProbes

    tmpSpks = orgSpks(animal,unit,expt,p,anaMode,dataFold);
    load(fullfile(physDir,animal,baseName,[baseName '_p' num2str(p) '_spkSort.mat']))


    probe(p).nSU = length(spkSort.unitinfo);
    probe(p).area = id.probes(p).area;
    for su = 1:probe(p).nSU

        spkTimes{p,su} = spkSort.spktimes(spkSort.unitid == su);

        for cond = 1:nConds
            for rep = 1:nReps

                trial = condTrials(rep,cond);
                spkTrain{p,su}(:,rep,cond) = ismember(nS(:,rep,cond),spkTimes{p,su}); 

                for bin = 1:nBins
                    binnedSpkTrain{p,su}(bin,rep,cond) = sum(spkTrain{p,su}((1:binSizeSam)+(binSizeSam*(bin-1)),rep,cond));
                end

%                 gausSpkTrain{p,su}(rep,:,cond) = conv(spkTrain{p,su}(rep,:,cond),kernel,'same');
%                 gausBinnedSpkTrain{p,su}(rep,:,cond) = conv(binnedSpkTrain{p,su}(rep,:,cond),kernel,'same');

            end
        end

    end

    if strcmp(areaA,probe(p).area)
        pA = p;
        %calc goodUnits
        [goodUnits{1},pVis{1}] = screenUnits(tmpSpks,anaMode,blank,'ranksum',0.01);
    elseif strcmp(areaB,probe(p).area)
        pB = p;
        %calc goodUnits
        [goodUnits{2},pVis{2}] = screenUnits(tmpSpks,anaMode,blank,'ranksum',0.01);
    end

end

%% Calculate cross corellograms 

for suA = 1:probe(pA).nSU
    for suB = 1:probe(pB).nSU
        disp([baseName ' ' probe(pA).area 'su#' num2str(suA) ' vs ' probe(pB).area 'su#' num2str(suB)])

        lambdaA = shiftdim(mean(sum(spkTrain{pA,suA},1)/timeLengthSec,2)); % mean firing rate (hz) responses to each condition, calculated as sum of spktrain matrix over time (rows) and mean over reps (col.s)
        lambdaB = shiftdim(mean(sum(spkTrain{pB,suB},1)/timeLengthSec,2));

        FR{1,suA} = shiftdim(sum(spkTrain{pA,suA},1)/timeLengthSec);
        FR{2,suB} = shiftdim(sum(spkTrain{pB,suB},1)/timeLengthSec);

        prefCondA{suA} = find(mean(sum(spkTrain{pA,suA},2),1) == max(mean(sum(spkTrain{pA,suA},2),1)));
        if length(prefCondA{suA}) > 1
            prefCondA{suA} = prefCondA{suA}(1);
        end
        prefCondB{suB} = find(mean(sum(spkTrain{pB,suB},2),1) == max(mean(sum(spkTrain{pB,suB},2),1)));
        if length(prefCondB{suB}) > 1
            prefCondB{suB} = prefCondB{suB}(1);
        end

% compute cross-corellogram ('ccg') for each unit pairing ({suA,suB}), for 
% each condition (col.), as the mean xcorr (rAB) across reps (col.s) 
        if xCorrMethod == 1
            
            for cond = 1:nConds
                for rep = 1:nReps
                    a = binnedSpkTrain{pA,suA}(:,rep,cond);
                    b = binnedSpkTrain{pB,suB}(:,rep,cond);
                    [rAB{suA,suB}(:,rep,cond),lag] = xcorr(a,b);
                end
                if normMethod == 1
                    ccgAB{suA,suB}(:,cond) = mean(rAB{suA,suB}(:,:,cond),2)./((nBins-abs(lag))*sqrt(mean(lambdaA)*mean(lambdaB)))';
                elseif normMethod == 2
                    ccgAB{suA,suB}(:,cond) = mean(rAB{suA,suB}(:,:,cond),2)./sqrt(mean(lambdaA)*mean(lambdaB));
                elseif normMethod == 3
                    ccgAB{suA,suB}(:,cond) = mean(rAB{suA,suB}(:,:,cond),2);
                end
            end
        
        elseif xCorrMethod == 2
    
%             for cond = 1:nConds
%                 for rep = 1:nReps
%     
%                     for tau = 1:length(lag)
%     
%                         a = binnedSpkTrain{pA,suA}(:,rep,cond);
%                         if lag(tau) == 0
%                             b = binnedSpkTrain{pB,suB}(:,rep,cond);
%                         elseif lag(tau)<0
%                             b = [binnedSpkTrain{pB,suB}((1+abs(lag(tau)):end),rep,cond);zeros(abs(lag(tau)),1)]; 
%                         elseif lag(tau)>0
%                             b = [zeros(abs(lag(tau)),1);binnedSpkTrain{pB,suB}((1:end-abs(lag(tau))),rep,cond)]; 
%                         end
%                         dpAB{suA,suB}(:,rep,cond,tau) = a.*b;
%     
%                     end
%     
%                     rAB{suA,suB}(:,rep,cond) = shiftdim(sum(dpAB{suA,suB}(:,rep,cond,:),1))./(sqrt(lambdaA(cond)*lambdaB(cond)))';
%                 end
%     
%                 if normMethod == 1
%                     ccgAB{suA,suB}(:,cond) = shiftdim(mean(sum(dpAB{suA,suB}(:,:,cond,:),1),2))./((nBins-abs(lag))*(sqrt(lambdaA(cond)*lambdaB(cond))))'; 
%                 elseif normMethod == 2
%                     ccgAB{suA,suB}(:,cond) = shiftdim(mean(sum(dpAB{suA,suB}(:,:,cond,:),1),2))./sqrt(lambdaA(cond)*lambdaB(cond));
%                 elseif normMethod == 3
%                     ccgAB{suA,suB}(:,cond) = shiftdim(mean(sum(dpAB{suA,suB}(:,:,cond,:),1),2));
%                 end
% 
% 
%             end
        
        end

    end
end




%% Calc CCG for neuron pairs in same area (area A vs area A; area B vs area B)


for suA1 = 1:probe(pA).nSU
    for suA2 = 1:probe(pA).nSU
        disp([baseName ' ' probe(pA).area 'su#' num2str(suA1) ' vs ' probe(pA).area 'su#' num2str(suA2)])

        lambdaA1 = shiftdim(mean(sum(spkTrain{pA,suA1},1)/timeLengthSec,2));
        lambdaA2 = shiftdim(mean(sum(spkTrain{pA,suA2},1)/timeLengthSec,2));

        if xCorrMethod == 1

            for cond = 1:nConds
                for rep = 1:nReps
                    a1 = binnedSpkTrain{pA,suA1}(:,rep,cond);
                    a2 = binnedSpkTrain{pA,suA2}(:,rep,cond);
                    [rAA{suA1,suA2}(:,rep,cond),lag] = xcorr(a1,a2);
                end

                if normMethod == 1
                    ccgAA{suA1,suA2}(:,cond) = mean(rAA{suA1,suA2}(:,:,cond),2)./((nBins-abs(lag))*sqrt(mean(lambdaA1)*mean(lambdaA2)))';
                elseif normMethod == 2
                    ccgAA{suA1,suA2}(:,cond) = mean(rAA{suA1,suA2}(:,:,cond),2)./sqrt(mean(lambdaA1)*mean(lambdaA2));
                elseif normMethod == 3
                    ccgAA{suA1,suA2}(:,cond) = mean(rAA{suA1,suA2}(:,:,cond),2);
                end

            end
        
        elseif xCorrMethod == 2
    
%             for cond = 1:nConds
%                 for rep = 1:nReps
%     
%                     for tau = 1:length(lag)
%     
%                         a1 = binnedSpkTrain{pA,suA1}(:,rep,cond);
%                         if lag(tau) == 0
%                             a2 = binnedSpkTrain{pA,suA2}(:,rep,cond);
%                         elseif lag(tau)<0
%                             a2 = [binnedSpkTrain{pA,suA2}((1+abs(lag(tau)):end),rep,cond);zeros(abs(lag(tau)),1)];
%                         elseif lag(tau)>0
%                             a2 = [zeros(abs(lag(tau)),1);binnedSpkTrain{pA,suA2}((1:end-abs(lag(tau))),rep,cond)];
%                         end
%                         dpAA(tau,:,rep,cond) = a1.*a2;
%     
%                     end
%     
%                 end
% 
%                 if normMethod == 1
%                     ccgAA{suA1,suA2}(:,cond) = shiftdim(mean(sum(dpAA(:,:,:,cond),2),3))./((nBins-abs(lag(tau)))*(sqrt(lambdaA1(cond).*lambdaA2(cond))))';
%                 elseif normMethod == 2
%                     ccgAA{suA1,suA2}(:,cond) = shiftdim(mean(sum(dpAA(:,:,:,cond),2),3))./(sqrt(lambdaA1(cond).*lambdaA2(cond)))';
%                 elseif normMethod == 3
%                     ccgAA{suA1,suA2}(:,cond) = shiftdim(mean(sum(dpAA(:,:,:,cond),2),3));
%                 end
%             end
        
        end

    end
end






for suB1 = 1:probe(pB).nSU
    for suB2 = 1:probe(pB).nSU
        disp([baseName ' ' probe(pB).area 'su#' num2str(suB1) ' vs ' probe(pB).area 'su#' num2str(suB2)])

        lambdaB1 = shiftdim(mean(sum(spkTrain{pB,suB1},1)/timeLengthSec,2));
        lambdaB2 = shiftdim(mean(sum(spkTrain{pB,suB2},1)/timeLengthSec,2));

        if xCorrMethod == 1

            for cond = 1:nConds
                for rep = 1:nReps
                    b1 = binnedSpkTrain{pB,suB1}(:,rep,cond);
                    b2 = binnedSpkTrain{pB,suB2}(:,rep,cond);
                    [rBB{suB1,suB2}(:,rep,cond),lag] = xcorr(b1,b2);
                end

                if normMethod == 1
                    ccgBB{suB1,suB2}(:,cond) = mean(rBB{suB1,suB2}(:,:,cond),2)./((nBins-abs(lag))*sqrt(mean(lambdaB1)*mean(lambdaB2)))';
                elseif normMethod == 2
                    ccgBB{suB1,suB2}(:,cond) = mean(rBB{suB1,suB2}(:,:,cond),2)./sqrt(mean(lambdaB1)*mean(lambdaB2));
                elseif normMethod == 3
                    ccgBB{suB1,suB2}(:,cond) = mean(rBB{suB1,suB2}(:,:,cond),2);    
                end
            end
        
        elseif xCorrMethod == 2
    
%             for cond = 1:nConds
%                 for rep = 1:nReps
%     
%                     for tau = 1:length(lag)
%     
%                         b1 = binnedSpkTrain{pB,suB1}(:,rep,cond);
%                         if lag(tau) == 0
%                             b2 = binnedSpkTrain{pB,suB2}(:,rep,cond);
%                         elseif lag(tau)<0
%                             b2 = [binnedSpkTrain{pB,suB2}((1+abs(lag(tau)):end),rep,cond);zeros(abs(lag(tau)),1)];
%                         elseif lag(tau)>0
%                             b2 = [zeros(abs(lag(tau)),1);binnedSpkTrain{pB,suB2}((1:end-abs(lag(tau))),rep,cond)];
%                         end
%                         dpBB(tau,:,rep,cond) = b1.*b2;
%     
%                     end
%     
%                 end
% 
%                 if normMethod == 1
%                     ccgBB{suB1,suB2}(:,cond) = shiftdim(mean(sum(dpBB(:,:,:,cond),2),3))./((nBins-abs(lag(tau)))*(sqrt(lambdaB1(cond).*lambdaB2(cond))))';
%                 elseif normMethod == 2
%                     ccgBB{suB1,suB2}(:,cond) = shiftdim(mean(sum(dpBB(:,:,:,cond),2),3))./(sqrt(lambdaB1(cond).*lambdaB2(cond)))';
%                 elseif normMethod == 3
%                     ccgBB{suB1,suB2}(:,cond) = shiftdim(mean(sum(dpBB(:,:,:,cond),2),3));
%                 end
% 
% 
%             end
        
        end

    end
end






%% Shuffle data

shuffs = [nchoosek(1:nReps,2);fliplr(nchoosek(1:nReps,2))];
[~,i] = sort(shuffs(:,1));
shuffs = shuffs(i,:);
for rep = 1:nReps %for each repeat, find non-matching repeat pairings

    shuffPairs(:,:,rep) = shuffs(shuffs(:,1)==rep,:); % matrix of non-matching repeat pairs organized by rep in the 3rd dim; size: nReps-1 x 2 x nReps

end


    for suA = 1:probe(pA).nSU % for each neuron pair...
        for suB = 1:probe(pB).nSU 
            disp([baseName ' (shuffle) ' probe(pA).area 'su#' num2str(suA) ' vs ' probe(pB).area 'su#' num2str(suB)])

            lambdaA = shiftdim(mean(sum(spkTrain{pA,suA},1)/timeLengthSec,2));
            lambdaB = shiftdim(mean(sum(spkTrain{pB,suB},1)/timeLengthSec,2));

            if xCorrMethod == 1

                for cond = 1:nConds % on a condition by condition basis...
                    for rep = 1:nReps % for each rep
                        for shuffRep = 1:nReps-1 % calculate the xcorr with all other reps of the same condition

                            sRep = shuffPairs(shuffRep,2,rep);
                            a = binnedSpkTrain{pA,suA}(:,rep,cond);
                            b = binnedSpkTrain{pB,suB}(:,sRep,cond);
                            [rShuff(:,shuffRep,rep),lag] = xcorr(a,b);

                        end
                    end

                    % shuffled ccg for each condition calculated as the mean xcorr output over
                    % non-matching pairs for each rep, then mean over reps,
                    % then normalized
                    if normMethod == 1
                        ccgShAB{suA,suB}(:,cond) = mean(mean(rShuff,2),3)./((nBins-abs(lag))*sqrt(mean(lambdaA)*mean(lambdaB)))';
                    elseif normMethod == 2
                        ccgShAB{suA,suB}(:,cond) = mean(mean(rShuff,2),3)./(sqrt(mean(lambdaA)*mean(lambdaB)));
                    elseif normMethod == 3
                        ccgShAB{suA,suB}(:,cond) = mean(mean(rShuff,2),3);
                    end
                end
            
            elseif xCorrMethod == 2
        
%                 for cond = 1:nConds
%                     for rep = 1:nReps
%     
%                         for tau = 1:length(lag)
%     
%                             a = shuffBinnedSpkTrain{pA,suA}(:,rep,cond);
%                             if lag(tau) == 0
%                                 b = shuffBinnedSpkTrain{pB,suB}(:,rep,cond);
%                             elseif lag(tau)<0
%                                 b = [shuffBinnedSpkTrain{pB,suB}((1+abs(lag(tau)):end),rep,cond);zeros(abs(lag(tau)),1)];
%                             elseif lag(tau)>0
%                                 b = [zeros(abs(lag(tau)),1);shuffBinnedSpkTrain{pB,suB}((1:end-abs(lag(tau))),rep,cond)];
%                             end
%                             dpShAB(tau,:,rep,cond,shuffNum) = a.*b;
%     
%                         end
%     
%                     end
% 
%                     if normMethod == 1
%                         ccgShAB{suA,suB}(:,cond,repeat) = shiftdim(mean(sum(dpShAB{suA,suB}(:,:,:,cond,repeat),2),3))./((nBins-abs(lag(tau)))*(sqrt(lambdaA(cond).*lambdaB(cond))))';
%                     elseif normMethod == 2
%                         ccgShAB{suA,suB}(:,cond,shuffNum) = shiftdim(mean(sum(dpShAB(:,:,:,cond,shuffNum),2),3))./(sqrt(lambdaA(cond).*lambdaB(cond)))';
%                     elseif normMethod == 3
%                         ccgShAB{suA,suB}(:,cond,repeat) = shiftdim(mean(sum(dpShAB{suA,suB}(:,:,:,cond,repeat),2),3));
%                     end
% 
% 
%                 end

            end

        end
    end


%% Calculate max(ccg)'s lag and diff. with shuff

for suA = 1:probe(pA).nSU

    for suB = 1:probe(pB).nSU

        t = lag.*binSizeSec;

        ccg = conv(mean(ccgAB{suA,suB},2,'omitnan'),kernel,'same');
        ccgSh = conv(mean(mean(ccgShAB{suA,suB},3,'omitnan'),2,'omitnan'),kernel,'same');
        if unique(isnan(ccg)) == 1 || unique(isnan(ccgSh)) == 1 
        
            peaksAB(suA,suB,2) = nan;
            peaksAB(suA,suB,1) = nan;
        
        else
            peak = max(ccg-ccgSh);
            peaksAB(suA,suB,2) = peak;
        
            tau = t(ccg-ccgSh==peaksAB(suA,suB,2));
            if length(tau)>1
                tau = min(tau);
            end
            peaksAB(suA,suB,1) = tau;
        end

        if prefCondA{suA}==trialInfo.blankId || prefCondB{suB}==trialInfo.blankId || length(trialInfo.dom)>1 || ~strcmp(trialInfo.dom,'ori')
            peaksAB(suA,suB,3) = nan;
        else
            peaksAB(suA,suB,3) =  abs(trialInfo.domval(prefCondA{suA})-trialInfo.domval(prefCondB{suB}));
        end

    end

end
 

%% Plot

if plt == 1

    for suA = 1:probe(pA).nSU
    for suB = 1:probe(pB).nSU
    
        figure;
        subplot(3,1,1)
            xp = [0 length(spkTrain{pA,suA}(:,1,1))]; xp = [xp fliplr(xp)]/sf+(stimTime-timeLengthSec);
            yp = [0 0 nReps nReps]+(nReps*(prefCondA{suA}-1));
            patch(xp,yp,'k','FaceColor','none','EdgeColor','r','LineWidth',1,'LineStyle','--');hold on
            xp = [0 length(spkTrain{pB,suB}(:,1,1))]; xp = [xp fliplr(xp)]/sf+(stimTime-timeLengthSec);
            yp = [0 0 nReps nReps]+(nReps*(prefCondB{suB}-1));
            patch(xp,yp,'k','FaceColor','none','EdgeColor','b','LineWidth',1,'LineStyle','--');hold on
            for cond = 1:nConds
                xp = [0 length(spkTrain{pA,suA}(:,1,1))]; xp = [xp fliplr(xp)];
                xp = xp/sf+(stimTime-timeLengthSec);
                yp = [0 0 nReps nReps]+(nReps*(cond-1));
                patch(xp,yp,'k','FaceColor',colors{cond},'EdgeColor','none','FaceAlpha',0.2);hold on
                [x,y] = find(spkTrain{pA,suA}(:,:,cond));
                y = y+(nReps*(cond-1));
                x = x/sf+(stimTime-timeLengthSec);
                plot(vertcat(x',x'),vertcat(y',y'-1),'r','LineWidth',2);hold on
                [x,y] = find(spkTrain{pB,suB}(:,:,cond));
                y = y+(nReps*(cond-1));
                x = x/sf+(stimTime-timeLengthSec);
                plot(vertcat(x',x'),vertcat(y',y'-1),'b','LineWidth',2)
                ticks(cond) =  (nReps/2)+(nReps*(cond-1));
                tickLabels{cond} = num2str(cond);
            end
            yticks(ticks)
            yticklabels(tickLabels)
            ylim([0,nConds*nReps+1])
            ylabel('stimulus condition')

% commenting this out since I am no longer making shuffled spkTrain matrix
%         subplot(6,1,2) 
%             for cond = 1:nConds
%                 xp = [0 length(spkTrain{pA,suA}(:,1,1))]; xp = [xp fliplr(xp)];
%                 xp = xp/sf+(stimTime-timeLengthSec);
%                 yp = [0 0 nReps nReps]+(nReps*(cond-1));
%                 patch(xp,yp,'k','FaceColor',colors{cond},'EdgeColor','none','FaceAlpha',0.2);hold on
%                 [x,y] = find(shuffSpkTrain{pA,suA}(:,:,cond));
%                 x = x/sf+(stimTime-timeLengthSec);
%                 y = y+(nReps*(cond-1));
%                 plot(vertcat(x',x'),vertcat(y',y'-1),'r','LineWidth',2);hold on
%                 [x,y] = find(shuffSpkTrain{pB,suB}(:,:,cond));
%                 y = y+(nReps*(cond-1));
%                 x = x/sf+(stimTime-timeLengthSec);
%                 plot(vertcat(x',x'),vertcat(y',y'-1),'b','LineWidth',2)
%                 ticks(cond) =  (nReps/2)+(nReps*(cond-1));
%                 tickLabels{cond} = num2str(cond);
%             end
%             yticks(ticks)
%             yticklabels(tickLabels)
%             ylim([0,nConds*nReps+1])
%             ylabel('stimulus condition (shuffled reps)')


        subplot(3,1,2);hold on
            for cond = 1:nConds
                if smooth == 1
                    plot(lag.*binSizeSec,conv(ccgAB{suA,suB}(:,cond),kernel,'same'),'Color',colors{cond});
                else
                    plot(lag.*binSizeSec,ccgAB{suA,suB}(:,cond),'Color',colors{cond});
                end
            end
        subplot(3,1,3);hold on
            if smooth == 1
                x = lag.*binSizeSec;
                y = conv(mean(ccgAB{suA,suB},2,'omitnan'),kernel,'same');
                ySh = conv(mean(ccgShAB{suA,suB},2,'omitnan'),kernel,'same');
                plot(x,y,'k')
                plot(x,ySh,'r--')
                plot(x,y-ySh,'g')
                plot(peaksAB(suA,suB,1),peaksAB(suA,suB,2),'bo')
                plot(x,repmat(peaksAB(suA,suB,2),1,length(x)),'b')
            else
                x = lag.*binSizeSec;
                y = mean(ccgAB{suA,suB},2,'omitnan');
                ySh = mean(ccgShAB{suA,suB},2,'omitnan');
                plot(x,y,'k')
                plot(x,ySh,'r--')
                plot(x,y-ySh,'g')
            end
            legend('mean ccg across conds','shuffled data ccg','shuff. corrected mean ccg','peak')
        sgtitle([probe(pA).area ' su#' num2str(suA) ' (red) vs ' probe(pB).area ' su#' num2str(suB) ' (blue)'])
        
    end
    end


    figure;hold on
    title(['ccgAA: ' probe(pA).area ' vs ' probe(pA).area])
    for suA1 = 1:probe(pA).nSU
    for suA2 = 1:probe(pA).nSU
    
            x = lag.*binSizeSec;
            y = conv(mean(ccgAA{suA1,suA2},2,'omitnan'),kernel,'same');
            plot(x,y,'k')
    
    end
    end
    

    figure;hold on
    title(['ccgBB: ' probe(pB).area ' vs ' probe(pB).area])
    for suB1 = 1:probe(pB).nSU
    for suB2 = 1:probe(pB).nSU
    
            x = lag.*binSizeSec;
            y = conv(mean(ccgBB{suB1,suB2},2,'omitnan'),kernel,'same');
            plot(x,y,'k')
    
    end
    end


    figure
    title('shuffle corrected ccg peaks')
    for suA = 1:probe(pA).nSU
        for suB = 1:probe(pB).nSU
            plot(peaksAB(suA,suB,1),peaksAB(suA,suB,2),'o');hold on
        end
    end
    xlabel('lag')
    ylabel('peak height (corr)')


end


end