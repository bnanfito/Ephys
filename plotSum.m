
function plotSum(dat)
plr = 1;
alignTC = 1;
svePlt = 0;

% dat = dat(dat.goodUnit,:);
nU = height(dat);
if svePlt == 0
    if nU>50
        uInd = randperm(nU,50);
    else
        uInd = 1:nU;
    end
elseif svePlt == 1
    uInd = 1:nU;
end

for u = uInd
    dom = dat.paramKey{u};
    isOri = (length(dom)==1 & sum(strcmp(dom,'ori'))==1) |...
            (length(dom)==2 & sum(strcmp(dom,'ori'))==1 & sum(contains(dom,'size'))==1 );
    isCon = (length(dom)==2 && sum(strcmp(dom,'ori'))==1 && sum(strcmp(dom,'contrast'))==1) ||...
            (length(dom)==3 && sum(strcmp(dom,'ori'))==1 && sum(strcmp(dom,'contrast'))==1 && sum(contains(dom,'size'))==1);

    figure;hold on
    nr = 2; nc = 2;
    predelay = 1;
    stimTime = 1;
    postdelay = 1;
    nTrials = max(dat.fr(u).trialNum(:));
    nConds = size(dat.cndKey{u},1);
    if strcmp(dat.area{u},'PSS')
        clr = 'r';
    elseif strcmp(dat.area{u},'V1')
        clr = 'b';
    else
        clr = 'k';
    end

    x = dat.spkTimes{u}(1,:);
    y = dat.spkTimes{u}(2,:);

    %PSTH
    subplot(nr,nc,1); hold on; box on
    binSize = 0.1;
    bins = -predelay:binSize:stimTime+postdelay;
    h = histogram(x,'BinEdges',bins);
    h.EdgeColor = 'none';
    h.FaceColor = clr;
    h.BinCounts = h.BinCounts/(nTrials*binSize);
    if ~isnan(dat.latency(u))
        xline(dat.latency(u),'--','LineWidth',2)
    end
    maxY = ceil(max(h.BinCounts)+(max(h.BinCounts)/10));
    patch([0 1 1 0],[0 0 maxY maxY],'k','EdgeColor','none','FaceAlpha',0.2)
    xlim([bins(1) bins(end)])
    ylim([0 maxY]);

    %RASTER (y = trial)
    subplot(nr,nc,3); hold on; box on
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

    %RASTER (y = condition)
    subplot(nr,nc,4); hold on; box on
    for t = 1:nTrials
        xT = x(y==t & idx);
        if isempty(xT)
            continue
        end
        yT = unique(dat.fr(u).trialCond(dat.fr(u).trialNum==t));
        plot(xT,yT,'k.')
    end
    xlim([-predelay stimTime+postdelay])
    if ismember('rBlank',dat.Properties.VariableNames)
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
    
    ttl = [dat.exptName{u} ' ' dat.area{u} ' ' dat.uInfo{u} ' ' num2str(dat.uID(u))];

    %TUNING CURVE
    if isOri
        ttl = [ttl  ' ldr=' num2str(dat.ldr(u))];
        x = dat.condition{u}(strcmp(dom,'ori'),:);
        y = dat.response{u};
        meanY = mean(y,'omitnan');
        x_pref = dat.oriPref(u);
        y_pref = dat.rPref(u);
        sem = std(y,'omitnan')/sqrt(size(y,1));
        if plr == 1
            subplot(nr,nc,2,polaraxes);hold on
    %         polarplot(deg2rad(x),y,'k.')
            polarplot(deg2rad([x x(1)]),mean([y y(:,1)],'omitnan') ,[clr '-'],'LineWidth',2)
            polarplot(deg2rad([x x(1)]),mean([y y(:,1)],'omitnan')+[sem sem(1)],[clr ':'],'LineWidth',1)
            polarplot(deg2rad([x x(1)]),mean([y y(:,1)],'omitnan')-[sem sem(1)],[clr ':'],'LineWidth',1)
    %         polarplot(repmat(deg2rad(x),2,1),mean(y,'omitnan')+([1;-1]*sem),clr,'LineWidth',2)
    %         polarplot(deg2rad(dat.oriPref(u)),dat.rPref(u),'r*')
    %         polarplot(deg2rad([0 dat.meanVec{u}.angDir]),[0 dat.meanVec{u}.magDir],clr,'LineWidth',2)
            polarplot(deg2rad([0 dat.meanVec{u}.angDir]),[0 dat.meanVec{u}.magDir*dat.ldr(u)],'g','LineWidth',2)
        else
            subplot(nr,nc,2); hold on; box on
            if alignTC == 1
                [x,meanY,i] = alignDirTuning(x,meanY);
                sem = sem(i);
                y = y(:,i);
            end
    %         plot(x,y,[clr '.'])
            plot(x,meanY,[clr '-'],'Marker','square','MarkerSize',7,'MarkerFaceColor',clr,'LineWidth',2)
            plot(repmat(x,2,1),meanY+([1;-1]*sem),clr,'LineWidth',2)
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
    elseif isCon
        ttl = [ttl  ' c50=' num2str(dat.c50(u))];
        x = dat.condition{u}(strcmp(dom,'contrast'),:,dat.oriPref(u));
        y = dat.response{u}(:,:,dat.oriPref(u));
        meanY = mean(y,'omitnan');
        sem = std(y,'omitnan')/sqrt(size(y,1));
        subplot(nr,nc,2); hold on; box on
        plot(x,meanY,[clr '-'],'Marker','square','MarkerSize',7,'MarkerFaceColor',clr,'LineWidth',2)
        plot(repmat(x,2,1),meanY+([1;-1]*sem),clr,'LineWidth',2)
    end
    
    sgtitle(ttl)
    figName = [dat.exptName{u} '_' dat.area{u} '_' dat.uInfo{u} '_' num2str(dat.uID(u))];
    if svePlt == 1
        figFold = fullfile('Y:\Brandon\data\Figures',dat.exptName{u}(1:5),dat.exptName{u});
        if ~isfolder(figFold)
            mkdir(figFold)
        end
        saveas(gcf,fullfile(figFold,figName),'fig')
%         saveas(gcf,fullfile(figFold,figName),'svg')
    end

    if svePlt == 1
        close gcf
    end


end

end