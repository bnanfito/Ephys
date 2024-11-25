

close all
plr = 1;

dat = data.pssaf;
dat = dat(dat.goodUnit,:);
nU = height(dat);
if nU>50
uInd = 1:floor(nU/50):nU;
else
uInd = 1:nU;
end
for u = uInd
    figure;hold on
    nr = 2; nc = 2;
    predelay = 1;
    stimTime = 1;
    postdelay = 1;
    nTrials = max(max(dat.fr(u).trialNum));
    nConds = max(max(dat.fr(u).trialCond));

    x = dat.spkTimes{u}(1,:);
    y = dat.spkTimes{u}(2,:);
    subplot(nr,nc,1);hold on
    binSize = 0.1;
    bins = -predelay:binSize:stimTime+postdelay;
    h = histogram(x,'BinEdges',bins);
    h.EdgeColor = 'none';
    h.BinCounts = h.BinCounts/(nTrials*binSize);
    if ~isnan(dat.latency(u))
        xline(dat.latency(u),'--','LineWidth',2)
    end
    xlim([bins(1) bins(end)])

    subplot(nr,nc,3);hold on
    idx = x>-1 & x<2;
    plot(x(idx),y(idx),'k.')
    trialExclude = dat.fr(u).trialNum(  isnan(dat.fr(u).bc)  );
    for t = trialExclude'
        patch([-predelay stimTime+postdelay stimTime+postdelay -predelay],[t-0.5 t-0.5 t+0.5 t+0.5],'r','EdgeColor','none','FaceAlpha',0.2)
    end
    if ~isnan(dat.latency(u))
        xline(dat.latency(u),'--','LineWidth',2)
    end
    xlim([-predelay stimTime+postdelay])
    ylim([0 nTrials+1])

    subplot(nr,nc,4);hold on
    for t = 1:nTrials
        xT = x(y==t & idx);
        if isempty(xT)
            continue
        end
        yT = unique(dat.fr(u).trialCond(dat.fr(u).trialNum==t));
        plot(xT,yT,'k.')
    end
    xlim([-predelay stimTime+postdelay])
    ylim([0 nConds+1])

    x = dat.condition{u}(strcmp(dat.paramKey{u},'ori'),:);
    y = dat.response{u};
    sem = std(y,'omitnan')/sqrt(size(y,1));
    if plr == 1
        subplot(nr,nc,2,polaraxes);hold on
%         polarplot(deg2rad(x),y,'k.')
        polarplot(deg2rad([x x(1)]),mean([y y(:,1)],'omitnan') ,'k-o')
        polarplot(repmat(deg2rad(x),2,1),mean(y,'omitnan')+([1;-1]*sem),'k')
        polarplot(deg2rad(dat.oriPref(u)),dat.rPref(u),'r*')
        polarplot(deg2rad([0 dat.meanVec{u}.angDir]),[0 dat.meanVec{u}.magDir],'k','LineWidth',2)
        polarplot(deg2rad([0 dat.meanVec{u}.angDir]),[0 dat.meanVec{u}.magDir*dat.ldr(u)],'g','LineWidth',2)
    else
        subplot(nr,nc,2);hold on
        plot(x,y,'k.')
        plot(x,mean(y,'omitnan'),'k-o')
        plot(repmat(x,2,1),mean(y,'omitnan')+([1;-1]*sem),'k')
    end


end