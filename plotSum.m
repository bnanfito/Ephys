
function plotSum(dat,plt,svePlt,figFold)
plr = 0;
alignTC = 0;

% dat = dat(dat.goodUnit,:);
nU = height(dat);
if plt == 1
    if nU>50
        uInd = randi(nU,1,50);
    else
        uInd = 1:nU;
    end
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
    nConds = size(dat.cndKey{u},1);

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
    patch([0 1 1 0],[0 0 max(y)+1 max(y)+1],'k','EdgeColor','none','FaceAlpha',0.2)
    xlabel('time (sec)')
    ylabel('trial #')

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
    if ~isempty(dat.rBlank{u})
        ylim([0 nConds+1])
    else
        ylim([0 nConds])
    end
    yticIdx = 1:2:nConds;
    yticks(1:2:nConds)
    yticklabels(num2str(dat.cndKey{u}(yticIdx,:)))
    patch([0 1 1 0],[0 0 max(y)+1 max(y)+1],'k','EdgeColor','none','FaceAlpha',0.2)
    xlabel('time (sec)')
    ylabel('condition')
    

    x = dat.condition{u}(strcmp(dat.paramKey{u},'ori'),:);
    y = dat.response{u};
    meanY = mean(y,'omitnan');
    x_pref = dat.oriPref(u);
    y_pref = dat.rPref(u);
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
        if alignTC == 1
            [x,meanY,i] = alignDirTuning(x,meanY);
            sem = sem(i);
            y = y(:,i);
        end
        plot(x,y,'k.')
        plot(x,meanY,'k-o')
        plot(repmat(x,2,1),meanY+([1;-1]*sem),'k')
        if alignTC == 1
            xlim([-180 180])
            xticks([-180 -90 0 90 180])
            xlabel('motion dir. relative to pref. (deg)')
        else
            plot(x_pref,y_pref,'r*')
            xlim([0 360])
            xticks([0 90 180 270])
            xlabel('motion direction (deg)')
        end
        ylabel('firing rate (Hz)')
    end
    
    
    ttl = [dat.exptName{u} ' ' dat.area{u} ' ' dat.uInfo{u} ' ' num2str(dat.uID(u))];
    sgtitle(ttl)
    figName = [dat.exptName{u} '_' dat.area{u} '_' dat.uInfo{u} '_' num2str(dat.uID(u))];
    if svePlt == 1
        if ~isfolder(figFold)
            mkdir(figFold)
        end
        saveas(gcf,fullfile(figFold,figName),'fig')
%         saveas(gcf,fullfile(figFold,figName),'svg')
    end

    if plt == 0
        close gcf
    end


end

end