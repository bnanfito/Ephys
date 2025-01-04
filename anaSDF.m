clear all
close all

animal = 'febn1';
unit = '000';
expt = '001';
probe = 2;
dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
anaMode = 'MU';

exptName = [animal '_u' unit '_' expt];
physDir = fullfile(dataFold,'Ephys');
exptDir = fullfile(physDir,animal,exptName);
load(fullfile(exptDir, [exptName '_id.mat'] ))
load(fullfile(exptDir, [exptName '_trialInfo.mat'] ))
load(fullfile(exptDir, [exptName '.analyzer']),'-mat')

varInfo = whos('probe');
if strcmp(varInfo.class,'char')
    area = probe;
    probe = find(strcmp({id.probes(:).area}',area));
end
clear varInfo
area = id.probes(probe).area;
sf = id.sampleFreq;
predelay = getparam('predelay',Analyzer);
stimTime = getparam('stim_time',Analyzer);
postdelay = getparam('postdelay',Analyzer);
trialL = predelay+stimTime+postdelay;
nTrials = length(trialInfo.triallist);
nConds = length(unique(trialInfo.triallist));
nDom = length(trialInfo.dom);
nReps = nTrials/nConds;
if isempty(trialInfo.blankId)
    blank = zeros(1,nConds)==1;
    blankTrial = [];
else
    blank = (1:nConds)==trialInfo.blankId;
    blankTrial = find(trialInfo.triallist==find(blank));
end
[sortTrialInd(:,1),sortTrialInd(:,2)] = sort(trialInfo.triallist);
nT = ((1:(trialL*sf))-(predelay*sf))/sf;
kernel = normpdf(-3:6/(sf*0.1):3);
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};

%% compute

[spks,trialExclude] = orgSpks(animal,unit,expt,probe,anaMode,dataFold);
nU = length(spks);

[~,lateSortIdx] = sort([spks.late]);
[~,zPosSortIdx] = sort([spks.zPos]);

s = cat(3,spks.train);
sdf = nan(size(s));
for u = 1:nU
   for t = 1:nTrials
       sT = s(:,t,u);
       sdf(:,t,u) = conv(sT,kernel,'same');
   end
end

figure;hold on
imagesc(squeeze(sum(sdf,2))')

% for u = 1:length(spks)
%     s = spks(u).stimCent(1,:);
%     g = normpdf(repmat(nT,length(s),1)',s,repmat(0.05,1,length(s)));
%     sdf(u,:) = mean(g,2,'omitnan');
%     sdf(u,:) = sdf(u,:)/max(sdf(u,:));
% end



% anaMode = 'MU';
% % load(['/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/training/Train_V1Cool/' anaMode '/ranksum & rPref above 2/Train_V1Cool_' anaMode 'dataSet.mat'])
% 
% tbl = data.v1bf;
% tbl = tbl(tbl.goodUnit,:);
% 
% [~,lateSortIdx] = sort(tbl.latency);
% for u = 1:height(tbl)
%     x = -1:0.01:2;
%     spks = tbl.spkTimes{u}(1,:);
%     g = normpdf(repmat(x,length(spks),1)',spks,repmat(0.05,1,length(spks)));
%     sdf(u,:) = mean(g,2,'omitnan');
%     sdf(u,:) = sdf(u,:)/max(sdf(u,:));
% end
% 
% figure; hold on
% imagesc(sdf(lateSortIdx,:))




