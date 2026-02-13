function [time,sdf,sdfC] = compute_sdf(data)

spkTimes = data.spkTimes{1};
trialKey = data.fr(1).trialNum;
cPref = find( data.condition{1}(strcmp(data.paramKey{1},'ori'),:) == data.oriPref(1) );
cNull = find( data.condition{1}(strcmp(data.paramKey{1},'ori'),:) == data.oriNull(1) );


nTrials = numel(trialKey);
nReps = size(trialKey,1);
nConds = size(trialKey,2);
repKey = repmat((1:nReps)',1,nConds);
condKey = repmat(1:nConds,nReps,1);

x = spkTimes(1,:);
y = spkTimes(2,:);
time = -2:0.001:4;
w = 0.05;
for t = 1:nTrials
    xT = x(y==t);
    for s = 1:length(xT)
        g{t}(:,s) = normpdf(time,xT(s),w);
    end
    sdf(:,t) = sum(g{t}/max(g{t},[],'all'),2)/(w);
end
for c = 1:nConds
    sdfC(:,c) = mean( sdf(:,ismember(1:nTrials,trialKey(:,c)))' ,'omitnan');
    sdfC_sem(:,c) = sem( sdf(:,ismember(1:nTrials,trialKey(:,c)))' );
end
prefC = 

SDF.nT = time;
SDF.mean = mean(sdf','omitnan');
SDF.sem = sem(sdf');
SDF.cMean = sdfC;
SDF.cSem = sdfC_sem;

figure;
subplot(2,2,1);hold on
plot(x,y,'.')
subplot(2,2,3);hold on
% plot(nT,[g{:}]/max([g{:}],[],'all'))
plot(time,SDF.mean,'k','LineWidth',2)
patch([SDF.nT fliplr(SDF.nT)],[SDF.mean+SDF.sem fliplr(SDF.mean-SDF.sem)],'k','EdgeColor','none','FaceAlpha',0.2)

clrs = hsv(nConds);

subplot(2,2,2);hold on
for t = 1:nTrials
    xT = x(y==t);
    if isempty(xT)
        continue
    end
    c = unique(condKey(trialKey==t));
    r = unique(repKey(trialKey==t));
    yT = c-((1/nReps)*(r-1));
    plot(xT,yT, 'k.')
end
for c = 1:nConds
    patch([min(time) max(time) max(time) min(time)],[0 0 1 1]+(c-1),clrs(c,:),'EdgeColor','none','FaceAlpha',0.2)
end
ylim([0 nConds])

subplot(2,2,4);hold on
for c = [cPref cNull]
    plot(SDF.nT,SDF.cMean(:,c),'Color',clrs(c,:),'LineWidth',1)
    patch([SDF.nT fliplr(SDF.nT)],[SDF.cMean(:,c)'+SDF.cSem(:,c)' fliplr(SDF.cMean(:,c)'-SDF.cSem(:,c)')],clrs(c,:),'EdgeColor','none','FaceAlpha',0.2)
end


end