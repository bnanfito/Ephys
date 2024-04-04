function  exptRaster(animal,unit,expt,anaMode,savePlot,dataFold)

% clear all
% close all
% animal = 'febh2';
% unit = '000';
% expt = '040';
% anaMode = 'MU';
% savePlot = 1;

baseName = [animal '_u' unit '_' expt];
disp(baseName)
physDir = fullfile(dataFold,'Ephys',animal,baseName);
figDir = fullfile(dataFold,'Figures',animal,baseName);

kernel = normpdf(-3:6/2000:3);
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
colors = repmat(colors,1,20);
load(fullfile(physDir,[baseName '_id.mat']),'id')
load(fullfile(physDir,[baseName '_trialInfo.mat']),'trialInfo')
if isfile(fullfile(physDir,[baseName '_temp.mat']))
    load(fullfile(physDir,[baseName '_temp.mat']))
end

sf = id.sampleFreq;
nTrials = length(trialInfo.triallist);
nEpochs = floor(length(trialInfo.eventTimes)/nTrials);
trialStart = downsample(trialInfo.eventTimes,nEpochs);
stimStart = downsample(trialInfo.eventTimes,nEpochs,1);
stimEnd = downsample(trialInfo.eventTimes,nEpochs,2);
trialEnd = downsample(trialInfo.eventTimes,nEpochs,3);

ampFile = fullfile(physDir, [baseName '_amplifier.dat']);
fileInfo = dir(ampFile);
% nSamps = fileInfo.bytes/(2*sum(vertcat(id.probes.nChannels)));
% nT = (1:nSamps)/sf;

sdf{length(id.probes)} = [];
for p = 1:length(id.probes)
%% Load Spks
    if strcmp(anaMode,'MU')
        load(fullfile(physDir,[baseName '_p' num2str(p) '_MUspkMerge.mat']))
        spks = MUspkMerge;
        MUThreshTrialData(fullfile(dataFold,'Ephys'),animal,unit,expt,p,'id',3,1,1)
        load(fullfile(physDir,[baseName '_p' num2str(p) '_MUThreshTrial.mat']),'MUThresh','MUThreshInfo')
        trialExclude = MUThreshInfo.trialExclude;
        clear MUspkMerge MUThresh MUThreshInfo
    elseif strcmp(anaMode,'SU')
        load(fullfile(physDir,[baseName '_p' num2str(p) '_spkSort.mat']))
        spks = spkSort;
        clear spkSort
    end

    if strcmp(id.probes(p).area,'V1')
        curColor = 'b';
    elseif strcmp(id.probes(p).area,'PSS')
        curColor = 'r';
    else
        continue
    end

    nSamps = max(spks.spktimes);
    nT = (1:nSamps)/sf;

%% Plot Probe
    fig{p} = figure('Position',[0 0 1800 1050]);
    
    if strcmp(anaMode,'MU')

        ax1 = subplot(2,1,1); hold on 
        for t = 1:nTrials 
            tStart = trialStart(t);
            sStart = stimStart(t);
            sEnd = stimEnd(t);
            tEnd = trialEnd(t);
            x = [sStart sEnd sEnd sStart]/sf;
            y = [0 0 65 65];
            if strcmp(anaMode,'MU') && ismember(t,find(trialExclude))
                patch([tStart tEnd tEnd tStart]/sf,y,'r','EdgeColor','none','FaceAlpha',0.2)
            end
            if ~isempty(trialInfo.blankId) && ismember(t, find(trialInfo.triallist==trialInfo.blankId)) 
                continue
            end
            patch(x,y,'k','EdgeColor','none','FaceAlpha',0.2) %%% stimulus %%%
        end
        x = spks.spktimes/sf;
        y = spks.detChSort;
        plot(x,y,[curColor '.']) %%% raster plot %%%
        xlabel('time (seconds)')
        ylim([0 max(y)+1])
        ylabel([id.probes(p).area ' det Ch sort'])

        ax2 = subplot(2,1,2); hold on
        hist = histogram(x,[0:0.1:nT(end)]);
        h = max(hist.Values)+10;
%         MUs = unique(spks.detCh);
%         sdf{p} = zeros(length(MUs),length(nT));
%         for u = 1:length(MUs)
%             disp(['computing sdf for MU#' num2str(u)])
%             sdf{p}(MUs(u),:) = conv(ismember(nT,x(spks.detCh==MUs(u))),kernel,'same');
%         end
%         h = max(mean(sdf{p}))+(max(mean(sdf{p}))/10);
        for t = 1:nTrials
            tStart = trialStart(t);
            sStart = stimStart(t);
            sEnd = stimEnd(t);
            tEnd = trialEnd(t);
            x = [sStart sEnd sEnd sStart]/sf;
            y = [0 0 h h];
            if strcmp(anaMode,'MU') && ismember(t,find(trialExclude))
                patch([tStart tEnd tEnd tStart]/sf,y,'r','EdgeColor','none','FaceAlpha',0.2)
            end
            if ~isempty(trialInfo.blankId) && ismember(t, find(trialInfo.triallist==trialInfo.blankId)) 
                continue
            end
            patch(x,y,'k','EdgeColor','none','FaceAlpha',0.2) %%% stimulus %%%
        end
%         plot(nT,mean(sdf{p}),curColor,'LineWidth',2); %%% SDF plot %%%
%         ylim([0 h])
        xticks(trialStart/sf);
        xticklabels([1:nTrials]);
        xlabel('trial')

        linkaxes([ax1 ax2],'x')

    elseif strcmp(anaMode,'SU')

%             uV = unique(v1Spks.unitid);
%             nSU = length(uV)-1;
%             for u = 1:nSU
%                 uID = uV(u);
%                 if uID == 0
%                     plot(v1Spks.spktimes(v1Spks.unitid==uID)/sf,v1Spks.detChSort(v1Spks.unitid == uID),'o','Color',[0.2 0.2 0.2],'MarkerSize',markerSz)
%                 else 
%                     plot(v1Spks.spktimes(v1Spks.unitid==uID)/sf,v1Spks.detChSort(v1Spks.unitid == uID),'o','Color',colors{u},'MarkerSize',markerSz)
%                 end
%             end

    end
    sgtitle(fig{p},[animal ' u' unit  ' e' expt ' p' num2str(p) ' (' id.probes(p).area ')'])

    if exist('temp','var') % plot temperature data if it exists
        yyaxis right
        stem(trialStart/sf,temp,'Color','green','LineWidth',2)
        ylim([-3 32])
        ylabel('temperature (C)')
    end

%% Save Plot
if savePlot==1
    if ~isdir(figDir)
        mkdir(figDir)
    end
    saveas(gcf,fullfile(figDir,[baseName '_p' num2str(p) '_exptRaster']),'jpeg')
%     for i = 1:3
%         if i==1
%             xlim([trialStart(1) trialEnd(5)]/sf)
%             fname = 'First5';
%         elseif i==2
%             xlim([trialStart(floor(nTrials/2)-5) trialEnd(floor(nTrials/2))]/sf)
%             fname = 'Mid5';
%         elseif i==3
%             xlim([trialStart(nTrials-5) trialEnd(nTrials)]/sf)
%             fname = 'Last5';
%         end
%         saveas(gcf,fullfile(figDir,[baseName '_p' num2str(p) '_exptRaster' fname]),'jpeg')
%     end
end

    

end


end