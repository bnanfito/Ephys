%plotMerge

function plotMerge(sumStats)

clrs = {'k','c','k'};
linStyls = {'-','-','--'};
plr = 1;
alignTC = 0;

nU = height(sumStats{1});
for u = 1:nU
    figure; hold on
    trialCount = 0;
    mergeId = [];
    for f = 1:length(sumStats)

        dat = sumStats{f};
        clr = clrs{f};
        linStyl = linStyls{f};
        mergeId = [mergeId dat.exptName{u}([8:10,12:14])];

        predelay = 1;
        stimTime = 1;
        postdelay = 1;
        nTrials = max(max(dat.fr(u).trialNum));
        nConds = size(dat.cndKey{u},1);
    
        x = dat.spkTimes{u}(1,:);
        y = dat.spkTimes{u}(2,:)+trialCount;
        if isempty(x)
            continue
        end

        %PSTH
        subplot(2,2,1);hold on; box on
        binSize = 0.1;
        bins = -predelay:binSize:stimTime+postdelay;
        h = histogram(x,'BinEdges',bins);
        h.EdgeColor = 'none';
        h.FaceColor = clr;
        h.BinCounts = h.BinCounts/(nTrials*binSize);
        if ~isnan(dat.latency(u))
            xline(dat.latency(u),'--','LineWidth',2,'Color',clr)
        end
        xlim([bins(1) bins(end)])


        %RASTER (trial)
        subplot(2,2,3);hold on; box on
        idx = x>-predelay & x<(stimTime+postdelay);
        plot(x(idx),y(idx),'.','Color',clr)
        trialExclude = dat.fr(u).trialNum(  isnan(dat.fr(u).bc)  );
        for t = trialExclude'
            patch([-predelay stimTime+postdelay stimTime+postdelay -predelay],[t-0.5 t-0.5 t+0.5 t+0.5],'r','EdgeColor','none','FaceAlpha',0.2)
        end
        if ~isnan(dat.latency(u))
            xline(dat.latency(u),'--','LineWidth',2,'Color',clr)
        end
        xlim([-predelay stimTime+postdelay])
        ylim([0 nTrials+trialCount])
        patch([0 1 1 0],[trialCount trialCount nTrials+trialCount nTrials+trialCount],'k','EdgeColor','none','FaceAlpha',0.2)
        xlabel('time (sec)')
        ylabel('trial #')


        %TUNING CURVE
        dom = dat.paramKey{u};
        isOri = (length(dom)==1 & sum(strcmp(dom,'ori'))==1) |...
                (length(dom)==2 & sum(strcmp(dom,'ori'))==1 & sum(contains(dom,'size'))==1 );
        isCon = (length(dom)==2 && sum(strcmp(dom,'ori'))==1 && sum(strcmp(dom,'contrast'))==1) ||...
                (length(dom)==3 && sum(strcmp(dom,'ori'))==1 && sum(strcmp(dom,'contrast'))==1 && sum(contains(dom,'size'))==1);
        if isOri
            y = dat.response{u};    
            meanY = mean(y,'omitnan');
            semY = std(y,'omitnan')/sqrt(size(y,1));
            x = dat.condition{u}(strcmp(dat.paramKey{u},'ori'),:);
            if plr == 1
%                 subplot(1,2,2,polaraxes);hold on;box on
% %                 polarplot(deg2rad(x),y,'.','Color',clr)
%                 polarplot(deg2rad([x x(1)]),mean([y y(:,1)],'omitnan') ,'-o','Color',clr)
%                 polarplot(repmat(deg2rad(x),2,1),mean(y,'omitnan')+([1;-1]*sem),'Color',clr)
%                 polarplot(deg2rad(dat.oriPref(u)),dat.rPref(u),'r*')
%                 polarplot(deg2rad([0 dat.meanVec{u}.angDir]),[0 dat.meanVec{u}.magDir],'k','LineWidth',2)
%                 polarplot(deg2rad([0 dat.meanVec{u}.angDir]),[0 dat.meanVec{u}.magDir*dat.ldr(u)],'g','LineWidth',2)
                
                subplot(1,2,2);hold on;box on; axis square
                %polar axes
                rtic = 0;
                while rtic<max(meanY)
                    rtic = rtic+5;
                    [axT,axR] = pol2cart(deg2rad(0:360),ones(1,361)*rtic);
                    plot(axT,axR,'k-')
                    rticMax = rtic;
                end
                ttic = 0;
                while ttic<360
                    [axT,axR] = pol2cart(deg2rad([ttic ttic]),[0 rticMax]);
                    plot(axT,axR,'k--')
                    ttic = ttic+45;
                end
                tmp = [meanY+semY fliplr(meanY-semY)];tmp(tmp<0) = 0;
                [errX,errY] = pol2cart(deg2rad([x fliplr(x)]),tmp);
                [rX,rY] = pol2cart(deg2rad(x),meanY);
                rX = [rX rX(1)];
                rY = [rY rY(1)];
                errX = [errX errX(1)];
                errY = [errY errY(1)];
                patch(errX,errY,clr,'EdgeColor','none','FaceAlpha',0.2)
                plot(errX,errY,clr)
                plot(rX,rY,clr,'LineWidth',2)

            else
                subplot(1,2,2);hold on;box on
                if alignTC == 1
                    [x,meanY,i] = alignDirTuning(x,meanY);
                    semY = semY(i);
                    y = y(:,i);
                else
                    x = [x 360];
                    y = [y y(:,1)];
                    meanY = [meanY meanY(1)];
                    semY = [semY semY(1)];
                end
%                 plot(x,y,'.','Color',clr)
                plot(x,meanY,'-square','Color',clr,'MarkerSize',7,'MarkerFaceColor',clr,'LineStyle',linStyl,'LineWidth',1)
                plot(repmat(x,2,1),meanY+([1;-1]*semY),'Color',clr,'LineWidth',1)
%                 plot(0:360,dat.gauss{u}.fit(0:360),'Color',clr,'LineWidth',2)
                if alignTC == 1
                    xlim([-180 180])
                    xticks([-180 -90 0 90 180])
                    xlabel('motion dir. relative to pref. (deg)')
                else
%                     plot(x_pref,y_pref,'r*')
                    xlim([0 360])
                    xticks([0 90 180 270])
                    xlabel('motion direction (deg)')
                end
                ylabel('firing rate (Hz)')
            end
        elseif isCon
            oriPref = dat.oriPref(u);
            y = dat.response{u}(:,:,oriPref);    
            meanY = mean(y,'omitnan');
            semY = std(y,'omitnan')/sqrt(size(y,1));
            x = dat.condition{u}(strcmp(dat.paramKey{u},'contrast'),:,oriPref);
            subplot(1,2,2);hold on;box on
            plot(x,meanY,'square','Color',clr,'MarkerSize',7,'MarkerFaceColor',clr)
            plot(repmat(x,2,1),meanY+([1;-1]*semY),'Color',clr,'LineWidth',2)
            plot([1:100],dat.nkFit{u}.fit(1:100),'Color',clr,'LineStyle',linStyl,'LineWidth',2)

        end

        trialCount = trialCount+nTrials;
    end
    ttl = [dat.exptName{u}(1:5) ' ' mergeId ' ' dat.area{u} ' ' dat.uInfo{u} ' ' num2str(dat.uID(u))];
    sgtitle(ttl)
end

end