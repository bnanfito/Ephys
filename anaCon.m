% anaCon
function [sumStats,spks] = anaCon(animal,unit,expt,probe,anaMode,dataFold,plt,svePlt)

%% Initialize

physDir = fullfile(dataFold,'Ephys');
figDir = fullfile(dataFold,'Figures');


%% Settings

exptName = [animal '_u' unit '_' expt];

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
kernel = normpdf(-3:6/2000:3);
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};

contInd = strcmp(trialInfo.dom,'contrast');
oriInd = strcmp(trialInfo.dom,'ori'); oris = unique(trialInfo.domval(:,oriInd));
sizeInd = contains(trialInfo.dom,'size');
posInd = strcmp(trialInfo.dom,'x_pos');
if nDom==2 && sum(contInd)==1 && sum(oriInd)==1 && sum(sizeInd)==0
    cndInclude = ~blank;
    cndInclude = cndInclude(1:size(trialInfo.domval,1));
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

dim{1} = 'response';
for d = 1:nDom
    dim{1+d} = trialInfo.dom{d};
end

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

    for o = 1:length(oris)
        curOri = trialInfo.domval(:,oriInd) == oris(o);
        if size(cndInclude,1) ~= size(curOri,1)
            curOri = curOri';
        end
        R{u}(:,:,o) = spks(u).fr.bc(:,cndInclude & curOri);
        C{u}(:,:,o) = trialInfo.domval(cndInclude & curOri,:)';
        maxR(u,o) = max(mean(R{u}(:,:,o),'omitnan'));
    end
    rPref(u) = max(maxR(u,:));
    oriPref(u) = find( maxR(u,:) == rPref(u) , 1 );
    rMean{u} = mean(R{u}(:,:,oriPref(u)),'omitnan');
    cont = C{u}(contInd,:,oriPref(u));

    [nk{u}.x,nk{u}.y,nk{u}.cF,nk{u}.resnorm,nk{u}.residuals,nk{u}.aic,nk{u}.bic] = nakaRush(rMean{u},cont,0,0);
    if nk{u}.cF>100
        cF(u) = nan;
    else
        cF(u) = nk{u}.cF;
    end

end

varNames = {'exptName','probe','area','uInfo','uID','spkTimes','latency','fr','paramKey','cndKey','oriPref','response','condition','rPref','nkFit','cF','rBlank'};
sumStats = table(exptID',probeID',areaID',{spks.info}',uID',{spks.stimCent}',vertcat(spks.late),vertcat(spks.fr),paramKey',cndKey',oriPref',R',C',rPref',nk',cF',Rblank','VariableNames',varNames);

[pVis] = visTest(sumStats);
sumStats.pVis = pVis;
goodUid = screenUnits(sumStats,anaMode);

if strcmp(anaMode,'MU')
    sumStats.xPos = vertcat(spks.xPos);
    sumStats.zPos = vertcat(spks.zPos);
end

%% Plot

if plt == 1

    for u = 1:height(sumStats)

        figure;

        for o = 1:length(oris)
            clr = colors{o};
   
            x = sumStats.spkTimes{u}(1,:);
            y = sumStats.spkTimes{u}(2,:);
            spkIdx = ismember(y,find(trialInclude));
            blankSpkIdx = ismember(y,blankTrial);
    
            subplot(2,2,1);hold on
            bins = -1:0.1:2;
            h = histogram(x(spkIdx),'BinEdges',bins);
            h.FaceColor = clr;
            h.EdgeColor = 'none';
            h = histogram(x(blankSpkIdx),'BinEdges',bins);
            h.FaceColor = 'k';
            h.EdgeColor = 'none';
            if ~isnan(sumStats.latency(u))
                xline(sumStats.latency(u),'--','LineWidth',2)
            end
            xlim([bins(1) bins(end)])
    
            subplot(2,2,3);hold on
            if isempty(x(spkIdx)) || isempty(y(spkIdx))
                text(0,0,'no spikes')
                continue
            end
            patch([0 1 1 0],[0 0 max(y(spkIdx))+1 max(y(spkIdx))+1],'k','EdgeColor','none','FaceAlpha',0.2)
            for t = 1:nTrials
                if ismember(t,find(trialExclude))
                    patch([-predelay stimTime+postdelay stimTime+postdelay -predelay],[t-0.5 t-0.5 t+0.5 t+0.5],'r','EdgeColor','none','FaceAlpha',0.2)
                end
            end
            plot(x(spkIdx),y(spkIdx),'.','Color',clr)
            plot(x(blankSpkIdx),y(blankSpkIdx),'k.')
            if ~isnan(sumStats.latency(u))
                xline(sumStats.latency(u),'--','LineWidth',2)
            end
            xlim([-1 2])
            ylim([0 max(y(spkIdx))+1])

        end

        xT = cont;
        yT = R{u}(:,:,oriPref(u));
        subplot(1,2,2);hold on
        plot(xT,rMean{u},'o','Color',clr)
        plot(repmat(xT,size(yT,1),1),yT,'.','Color',clr)
        plot(nk{u}.x,nk{u}.y)
        xline(nk{u}.cF,'g--')
        yline(0,'k')

        ttl = [sumStats.uInfo{u} '#' num2str(sumStats.uID(u))];
        if ~ismember(u,find(goodUid))
            ttl = ['(BAD) ' ttl];
        end
        sgtitle(ttl)
        clear ttl

    end

end

%     figure;
%     cdfplot(sumStats.cF(sumStats.goodUnit))







end
