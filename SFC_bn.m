%Spike Field Coherence (SFC) analysis
%written by Brandon Nanfito
clear all
close all

%% Settings

if ispc
%     dataFold = 'D:\data'; 
%     dataFold = 'C:\Users\brand\Documents\data';
    dataFold = 'F:\Brandon\data';
elseif ismac
%     dataFold = '/Volumes/Lab drive/Brandon/data';
%     dataFold = '/Volumes/Lab drive/Brandon/VSS2024/data';
    dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end

saveFig=0;

animal = 'febl0'; 
unitID = '000'; 
expt = '005';
probe = 1;
unit = 6;
exptName =[animal '_u' unitID '_' expt];
exptPath = fullfile(dataFold,'Ephys',animal,exptName);


cd(exptPath);
load(fullfile(exptPath,[exptName '.analyzer']),'-mat')
load(fullfile(exptPath,[exptName '_id.mat']))
load(fullfile(exptPath,[exptName '_p' num2str(probe) '_spkSort.mat']))

sf=id.sampleFreq;
nProbes = length(id.probes);
nCh=sum(vertcat(id.probes.nChannels));
nUnits = length(spkSort.unitinfo);


%% compute channel neighbor matrix

chXpos{nProbes} = [];
for p = 1:nProbes
    if strcmp(id.probes(p).type,'64F')
        probe_64F;close all
        chXpos{p}=s.x;chZpos{p}=s.z;
    elseif strcmp(id.probes(p).type,'64D')
        probe_64D;close all
        chXpos{p}=s.x;chZpos{p}=s.z;
    end
end

dist = NaN(64,64,nProbes);
neigh = NaN(64,64,nProbes);
for p = 1:nProbes
    for ch = 1:64
        for di = 1:64
            dist(ch,di,p) = sqrt(((chXpos{p}(ch)-chXpos{p}(di))^2)+((chZpos{p}(ch)-chZpos{p}(di))^2));
            [D,N] = sort(dist(ch,:,p));
            neigh(ch,:,p) = N;
%             if ch == detCh && p == 1
%                 chInd = N;
%             end
        end
    end
end

%% extract STA for each unit
hp = 1.5;
lp = 150;

figure;hold on
for u = 1:nUnits
    if u == 6
        continue
    end
    disp(['unit: ' num2str(u)])
    spkTimes{u} = spkSort.spktimes(spkSort.unitid == u);

    for spk = 1:length(spkTimes{u})
        disp(['spk #' num2str(spk) '/' num2str(length(spkTimes{u}))])
        spkT = spkTimes{u}(spk);
        [~,dataFilt] = readRawData(animal,unitID,expt,spkT-(sf*1),spkT+(sf*1),[hp lp],dataFold);
        if spk>1 & size(dataFilt,2)~=size(lfp{u},2)
            continue
        end
        
        if probe == 1
            lfp{u}(spk,:) = mean(dataFilt( id.probes(probe).nChannels+1:end ,:));
        elseif probe == 2
            lfp{u}(spk,:) = mean(dataFilt( 1:id.probes(probe-1).nChannels ,:));
        end
        fclose all;
    end
    STA(u,:) = mean(lfp{u});
    plot(STA(u,:),'LineWidth',2)

end