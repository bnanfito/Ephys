% anaOri
% function [sumStats] = anaSF(animal,unit,expt,probe,anaMode,dataFold,plt,svePlt)

clear all
% close all
animal = 'febo5';
unit = '000';
expt = '008';
probe = 1;
anaMode = 'MU';
dataFold = '/Volumes/Lab drive/Brandon/data';


%% Initialize

physDir = fullfile(dataFold,'Ephys');
figDir = fullfile(dataFold,'Figures');


%% Settings

exptName = [animal '_u' unit '_' expt];

plr = 1;
alignBit = 0;
visTest = 'ranksum';
alpha = 0.01;


%% Load Data

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
sampFreq = id.sampleFreq;
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
nT = ((1:(trialL*sampFreq))-(predelay*sampFreq))/sampFreq;
kernel = normpdf(-3:6/2000:3);
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};

sfInd = strcmp(trialInfo.dom,'s_freq');
oriInd = strcmp(trialInfo.dom,'ori');
sizeInd = contains(trialInfo.dom,'size');
posInd = strcmp(trialInfo.dom,'x_pos');
if nDom==2 && sum(sfInd)==1 && sum(oriInd)==1
    cndInclude = ~blank;
elseif nDom == 3 && sum(sizeInd)==1 
    sizes = unique(trialInfo.domval(:,sizeInd));
    ffIdx = sizes>100;
    hemiIdx =  sizes<100;
    cndInclude = trialInfo.domval(:,sizeInd) == sizes(hemiIdx);
elseif nDom == 3 && sum(posInd)==1 
    pos = unique(trialInfo.domval(:,posInd));
    ffIdx = pos==min(pos);
    hemiIdx = pos==max(pos);
    cndInclude = trialInfo.domval(:,sizeInd) == pos(hemiIdx);
end
trialInclude = ismember(trialInfo.triallist,find(cndInclude));


%% Compute

[spks,trialExclude] = orgSpks(animal,unit,expt,probe,anaMode,dataFold);
nU = length(spks);

for u = 1:nU % u indexes a unit (column) in structure spks

    exptID{u} = exptName;
    uID(u) = spks(u).unitId;
    probeID(u) = probe;
    areaID{u} = area;
    paramKey{u} = trialInfo.dom;
    cndKey{u} = trialInfo.domval;

    if ~isempty(trialInfo.blankId)
        Rblank{u} = spks(u).fr.stim(:,blank);
    end

    c = trialInfo.domval(cndInclude,:)';


    dir = unique(c(oriInd,:));
    sf = unique(c(sfInd,:));

    for d = 1:length(dir)

        R{u}(:,:,d) = spks(u).fr.bc(:, c(oriInd,:)==dir(d));
        C{u}(:,:,d) = c(:,c(oriInd,:)==dir(d));

        rPref(u,d) = max(mean(R{u}(:,:,d),'omitnan'));
        sfPref(u,d) = C{u}(sfInd, find(mean(R{u}(:,:,d),'omitnan') == rPref(u,d),1) , d);

    end

    opIdx = rPref(u,:) == max(rPref(u,:));
    if sum(opIdx) == 1
        oriPref(u) = dir( opIdx );
    else 
        dirMean = squeeze(mean(mean(R{u}(:,:,:),1,'omitnan'),2,'omitnan'));
        if length(unique(dirMean))==1
            oriPref(u) = nan;
        else
            oriPref(u) = dir( dirMean == max( dirMean ) );
        end
    end


end

varNames = {'exptName','probe','area','uInfo','uID','spkTimes','latency','fr','paramKey','cndKey','response','condition','rPref','sfPref','oriPref'};
sumStats = table(exptID',probeID',areaID',{spks.info}',uID',{spks.stimCent}',vertcat(spks.late),vertcat(spks.fr),paramKey',cndKey',R',C',rPref,sfPref,oriPref','VariableNames',varNames);

if exist('Rblank','var')
    sumStats.rBlank = Rblank';
end

[goodUnit,pVis] = screenUnits(sumStats,anaMode,blank,visTest,alpha);
sumStats.goodUnit = goodUnit';
sumStats.pVis = pVis;


%% plotting

for u = 1:nU

    if ~sumStats.goodUnit(u)
        continue
    end

    opIdx = squeeze(mean(sumStats.condition{u}(oriInd,:,:),2)) == sumStats.oriPref(u);
%     figure;
    y(u,:) = mean( sumStats.response{u}(:,:,opIdx) ,1,'omitnan');
    y(u,:) = y(u,:)./max(y(u,:));
    x = sumStats.condition{u}(sfInd,:,opIdx);
%     plot(x,y);



end

figure;
plot(x,mean(y))



% end