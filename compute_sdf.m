function [SDF] = compute_sdf(data,varargin)

spkTimes = data.spkTimes{1};
conds = data.condition{1}(strcmp(data.paramKey{1},'ori'),:);
cPref = find( conds == data.oriPref(1) );
cNull = find( conds == data.oriNull(1) );
cOrth1 = find( conds == mod(conds(cPref)+90,360) );
cOrth2 = find( conds == mod(conds(cPref)-90,360) );
nConds = length(conds);
if length(data.paramKey{1})>1
    rA = data.response{1}; rA = rA(~isnan(rA));
    rB = data.fr(1).bc(:,1:nConds); rB= rB(~isnan(rB));
    rC = data.fr(1).bc(:,(1:nConds)+nConds); rC = rC(~isnan(rC));
    if unique(rB == rA)==1
        cndsInclude = 1:nConds;
    elseif unique(rC == rA)==1
        cndsInclude = (1:nConds)+nConds;
    end
else
    cndsInclude = 1:nConds;
end
trialKey = data.fr(1).trialNum(1:5,cndsInclude);
trialInclude = sort(trialKey(:));
nTrials = numel(trialKey);
nReps = size(trialKey,1);
repKey = repmat((1:nReps)',1,nConds);
condKey = repmat(1:nConds,nReps,1);

spkT = spkTimes(1,:);
trialId = spkTimes(2,:);
dt = 0.001;
time = -1:dt:2;
w = 0.020;
for t = 1:nTrials
    tIdx = trialInclude(t);
    if ismember(tIdx,trialId)
        xT = spkT(trialId==tIdx);
        for s = 1:length(xT)
            g{t}(:,s) = normpdf(time,xT(s),w)*dt;
        end
%         sdf(:,tId) = sum(g{tId}/max(g{tId},[],'all'),2)/(w);
        sdf(:,t) = sum(g{t},2)/dt;
    else
        sdf(:,t) = zeros(size(time));
    end
end
for c = cndsInclude
    sdf_cMean(c,:) = mean( sdf(:,ismember(1:nTrials,trialKey(:,c)))' ,'omitnan');
    sdf_cSem(c,:) = sem( sdf(:,ismember(1:nTrials,trialKey(:,c)))' );
end

SDF.time = time;
SDF.trial = sdf';
SDF.mean = mean(SDF.trial,'omitnan');
SDF.sem = sem(SDF.trial);
SDF.ci95 = confInt(SDF.trial);
SDF.cMean = sdf_cMean;
SDF.cSem = sdf_cSem;
SDF.cPref = cPref;

%% Plot

if ~isempty(varargin) && varargin{1}==1

    figure;
    tiledlayout(2,2)
    clrs = hsv(nConds);
    
    ax(1) = nexttile;hold on
    plot(spkT,trialId,'.')
    xlim([min(SDF.time) max(SDF.time)])
    ylim([0 nTrials])
    
    ax(2) = nexttile;hold on
    for c = 1:nConds
        patch([0 1 1 0],[0 0 1 1]+(c-1),clrs(c,:),'EdgeColor','none','FaceAlpha',0.2)
    end
    for t = 1:nTrials
        tIdx = trialInclude(t);
        xT = spkT(trialId==tIdx);
        if isempty(xT)
            continue
        end
        c = unique(condKey(trialKey==tIdx));
        r = unique(repKey(trialKey==tIdx));
        yT = c-((1/nReps)*(r-1));
        plot(xT,yT, 'k.')
    end
    xlim([min(SDF.time) max(SDF.time)])
    ylim([0 nConds])
    
    ax(3) = nexttile;hold on
    % plot(nT,[g{:}]/max([g{:}],[],'all'))
    plot(SDF.time,SDF.mean,'k','LineWidth',2)
    patch([SDF.time fliplr(SDF.time)],[SDF.mean+SDF.sem fliplr(SDF.mean-SDF.sem)],'k','EdgeColor','none','FaceAlpha',0.2)
    
    ax(4) = nexttile;hold on
    for c = [cPref,cNull]
        plot(SDF.time,SDF.cMean(c,:),'Color',clrs(c,:),'LineWidth',1)
        patch([SDF.time fliplr(SDF.time)],[SDF.cMean(c,:)+SDF.cSem(c,:) fliplr(SDF.cMean(c,:)-SDF.cSem(c,:))],clrs(c,:),'EdgeColor','none','FaceAlpha',0.2)
    end
    linkaxes(ax,'x')
end

end