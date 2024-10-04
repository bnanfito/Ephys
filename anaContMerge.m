function anaContMerge(animal,unit,expt,probe,anaMode,dataFold,plt,svePlt)

    close all
    
    colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};


    for f = 1:length(expt)

        curUnit = unit{f};
        curExpt = expt{f};
        [sumStats{f}] = anaCon(animal,curUnit,curExpt,probe,anaMode,dataFold,0,svePlt);

    end

    nU = height(sumStats{1});
    for u = 1:nU

        figure; hold on
        
        trialCount = 0;
        for f = 1:length(expt)

            clr = colors{f};

            fileBase = [animal '_u' curUnit '_' curExpt];
            load(fullfile(dataFold,'Ephys',animal,fileBase,[fileBase '_trialInfo.mat']),'trialInfo')
            nTrials = length(trialInfo.triallist);
            
            subplot(2,2,1); hold on
            x = sumStats{f}.spkTimes{u}(1,:);
            y = sumStats{f}.spkTimes{u}(2,:); y = y+trialCount;
            trialCount = trialCount+nTrials;
            blankSpkIdx = ismember(y,trialInfo.blankId);
    
            subplot(2,2,1);hold on
            bins = -1:0.1:2;
            h = histogram(x,'BinEdges',bins);
            h.FaceColor = clr;
            h.EdgeColor = 'none';
            h = histogram(x(blankSpkIdx),'BinEdges',bins);
            h.FaceColor = 'k';
            h.EdgeColor = 'none';
            if ~isnan(sumStats{f}.latency(u))
                xline(sumStats{f}.latency(u),'--','LineWidth',2)
            end
            xlim([bins(1) bins(end)])
    
            subplot(2,2,3);hold on
            if isempty(x) || isempty(y)
                text(0,0,'no spikes')
                continue
            end
            patch([0 1 1 0],[0 0 max(y)+1 max(y)+1],'k','EdgeColor','none','FaceAlpha',0.2)
            plot(x,y,'.','Color',clr)
            plot(x(blankSpkIdx),y(blankSpkIdx),'k.')
            if ~isnan(sumStats{f}.latency(u))
                xline(sumStats{f}.latency(u),'--','LineWidth',2)
            end
            xlim([-1 2])
            ylim([0 max(y)+1])

            subplot(1,2,2);hold on
            xT = sumStats{f}.condition{u}(strcmp(sumStats{f}.paramKey{u},'contrast'),:,sumStats{f}.oriPref(u));
            yT = sumStats{f}.response{u}(:,:,sumStats{f}.oriPref(u));
            rMean = mean(yT,'omitnan');
            plot(xT,rMean,'o','Color',clr)
            plot(repmat(xT,size(yT,1),1),yT,'.','Color',clr)
            plot(sumStats{f}.nkFit{u}.x,sumStats{f}.nkFit{u}.y,'Color',clr)
            xline(sumStats{f}.nkFit{u}.cF,'g--')
            yline(0,'k')

        end

        ttl = [sumStats{f}.uInfo{u} '#' num2str(sumStats{f}.uID(u))];
        if ~ismember(u,find(sumStats{f}.goodUnit))
            ttl = ['(BAD) ' ttl];
        end
        sgtitle(ttl)
        clear ttl

    end


end
