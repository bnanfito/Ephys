% anaOri
function [sumStats,spks] = anaOri(animal,unit,expt,probe,anaMode,dataFold,plt,svePlt,varargin)

%% Initialize

physDir = fullfile(dataFold,'Ephys');
figDir = fullfile(dataFold,'Figures');


%% Settings

exptName = [animal '_u' unit '_' expt];

plr = 1;
alignBit = 0;

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

if isempty(varargin)
    probeId = probe;
else
    % varargin{1} should be a merge file id (1,2,3..)
    probeId = [num2str(probe) '_' num2str(varargin{1})];
end

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
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],...
    [0.6350 0.0780 0.1840]};

oriInd = strcmp(trialInfo.dom,'ori');
sizeInd = contains(trialInfo.dom,'size');
posInd = strcmp(trialInfo.dom,'x_pos');
if nDom==1 && sum(oriInd)==1 && sum(sizeInd)==0
    cndInclude = ~blank;
elseif nDom == 2 && sum(sizeInd)==1 
    sizes = unique(trialInfo.domval(:,sizeInd));
    ffIdx = sizes>100;
    hemiIdx =  sizes<100;
    cndInclude = trialInfo.domval(:,sizeInd) == sizes(hemiIdx);
elseif nDom == 2 && sum(posInd)==1 
    pos = unique(trialInfo.domval(:,posInd));
    ffIdx = pos==min(pos);
    hemiIdx = pos==max(pos);
    cndInclude = trialInfo.domval(:,sizeInd) == pos(hemiIdx);
end
trialInclude = ismember(trialInfo.triallist,find(cndInclude));

%% Compute 

[spks,trialExclude] = orgSpks(animal,unit,expt,probeId,anaMode,dataFold);
nU = length(spks);

for u = 1:nU % u indexes a unit (column) in structure spks

    exptID{u} = exptName;
    uID(u) = spks(u).unitId;
    probeID(u) = probe;
    areaID{u} = area;
    paramKey{u} = trialInfo.dom;
    cndKey{u} = trialInfo.domval;

    if ~isempty(trialInfo.blankId)
        rBlank{u} = spks(u).fr.stim(:,blank);
    end

    R{u}(:,:) = spks(u).fr.bc(:,cndInclude);
    C{u}(:,:) = trialInfo.domval(cndInclude,:)';

    dir = C{u}(oriInd,:);

    if length(dir)>=4
        %compute metrics in orientation space
        ori = mod(dir,180);
        oris = unique(ori);
        for o = 1:length(oris)
            idx = ori == oris(o);
            rMean_ori(o) = mean(R{u}(:,idx),'all','omitnan');
        end
        rPref_ori = max(rMean_ori);
        if sum(rMean_ori==rPref_ori)>1
            rIn = rMean_ori;
            pks = rIn == rPref_ori;
            rConv = rIn([end 1:end 1]);
            rConv = conv(rConv,ones(1,3)*(1/3),'same');
            rConv = rConv(1+1:end-1);
            rConv(~pks) = 0;
            cPref_ori = oris(find(rConv==max(rConv),1,'first'));
            clear rIn pks rConv 
        else
            cPref_ori = oris(rMean_ori==rPref_ori);
        end
        cNull_ori = mod(cPref_ori+90,180);
        rNull_ori = rMean_ori(oris==cNull_ori);
        osi(u,1) = abs(rPref_ori-rNull_ori)/rPref_ori;
        clear rMean_ori rPref_ori rNull_ori

        %compute metrics in direction space
        rMean = mean(R{u}(:,:),'omitnan');
        rPref(u) = max(rMean);
        if sum(rMean==rPref(u))>1
            rIn = rMean;
            pks = rIn == rPref(u);
            rConv = rIn([end 1:end 1]);
            rConv = conv(rConv,ones(1,3)*(1/3),'same');
            rConv = rConv(1+1:end-1);
            rConv(~pks) = 0;
            cPref(u) = dir(find(rConv==max(rConv),1,'first'));
            clear rIn pks rConv 
        else
            cPref(u) = dir(rMean==rPref(u));
        end
        cNull(u) = mod(cPref(u)+180,360);
        rNull(u) = rMean(dir==cNull(u));
        dsi(u,1) = abs(rPref(u)-rNull(u))/rPref(u);
    
        %compute resultant vector
        mv{u} = meanvec(dir,rMean);
        ldr(u,1) = mv{u}.ldr;
        lor(u,1) = mv{u}.lor;
    
%         %double gaussian fit
%         G{u} = dirGauss(rMean,dir,0);
    elseif length(dir)==2
        %compute direction preference index
        rMean = mean(R{u},'omitnan');
        rPref(u) = max(rMean);
        if sum(rMean==rPref(u))>1
            cPref(u) = dir(find(rMean==rPref(u),1,'first'));
        else
            cPref(u) = dir(rMean==rPref(u));
        end
        cNull(u) = mod(cPref(u)+180,360);
        rNull(u) = rMean(dir==cNull(u));
        dpi(u,1) = abs(diff(rMean))/rPref(u);
    end

end

varNames = {'exptName','probe','area','uInfo','uID','spkTimes','latency','fr','paramKey','cndKey','response','condition','rPref','oriPref','rNull','oriNull'};
sumStats = table(exptID',probeID',areaID',{spks.info}',uID',{spks.stimCent}',vertcat(spks.late),vertcat(spks.fr),paramKey',cndKey',R',C',rPref',cPref',rNull',cNull','VariableNames',varNames);

if exist('rBlank','var')
    sumStats.rBlank = rBlank';
end

if exist('dpi','var')
    sumStats.dpi = dpi;
else
    sumStats.meanVec = mv';
    sumStats.dsi = dsi;
    sumStats.ldr = ldr;
    sumStats.osi = osi;
    sumStats.lor = lor;
end


[pVis] = visTest(sumStats);
sumStats.pVis = pVis;
goodUnit = screenUnits(sumStats,anaMode);

if strcmp(anaMode,'MU')
    sumStats.xPos = vertcat(spks.xPos);
    sumStats.zPos = vertcat(spks.zPos);
end


%% Plot

if plt == 1

    for u = 1:height(sumStats)

        if ~ismember(u,find(goodUnit))
            continue
        end
        
        figure;
        clr = colors{1};

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



        xT = sumStats.condition{u}(oriInd,:);
        yT = sumStats.response{u}(:,:);
        yMean = mean(yT,'omitnan');
        if alignBit == 1
            [xT,yMean,i] = alignDirTuning(xT,yMean);
            yT{u} = yT{u}(:,i);
        elseif exist('dsi','var')
            xT = [xT xT(1)+360];
            yT = [yT yT(:,1)];
            yMean = [yMean yMean(1)];
        end
        sem = std(yT,'omitnan')/sqrt(size(yT,1));

        if plr == 1
            subplot(1,2,2,polaraxes);hold on
            xT = deg2rad(xT);
            polarplot(xT,yMean,'-','Color',clr);
            hold on;
%                 polarplot(xP,yP,'.','Color',clr)
            polarplot(repmat(xT,2,1),yMean+([1;-1]*sem),'Color',clr,'LineWidth',2)
            if exist('dsi','var')
                polarplot(repmat(deg2rad(sumStats.meanVec{u}.angDir),2,1),...
                    [0 sumStats.meanVec{u}.magDir],'k','LineWidth',2)
                polarplot(repmat(deg2rad(sumStats.meanVec{u}.angDir),2,1),...
                    [0 sumStats.meanVec{u}.ldr*sumStats.meanVec{u}.magDir],'g','LineWidth',2)
            end
            polarplot(deg2rad(sumStats.oriPref(u)),sumStats.rPref(u),'r*')
        else
            subplot(1,2,2);hold on
            plot(xT,yT','.','Color',clr)
            plot(repmat(xT,2,1),yMean+([1;-1]*sem),'Color',clr)
            plot(xT,yMean,'o','Color',clr);
            if exist('dsi','var')
                plot(repmat(sumStats.meanVec{u}.angDir,2,1),[0;sumStats.meanVec{u}.magDir],'k','LineWidth',2)
                plot(repmat(sumStats.meanVec{u}.angDir,2,1),[0;sumStats.meanVec{u}.ldr*sumStats.meanVec{u}.magDir],'g','LineWidth',2)
            end
            plot(sumStats.oriPref(u),sumStats.rPref(u),'r*')
        end

        ttl = [exptName ' p' num2str(probe) ' (' area ') ' anaMode '#' num2str(uID(u))];
        if ~ismember(u,find(goodUnit))
            ttl = [ttl '(BAD UNIT)'];
        end
        sgtitle(ttl)
        if svePlt == 1 && ismember(u,find(goodUnit))
            figFileDir = fullfile(figDir,animal,exptName);
            if ~isfolder(figFileDir)
                mkdir(figFileDir)
            end
            figFileName = fullfile(figFileDir,[exptName '_p' num2str(probe) '_' anaMode num2str(uID(u))]);
            saveas(gcf,figFileName)
        end
    end

    
end



end
