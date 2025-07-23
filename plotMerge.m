%plotMerge

function [sumStats] = plotMerge(animalId, mergeId, probe, dataFold, plt)

anaMode = 'SU';

mergeName = [animalId '_uMMM_' mergeId];
splitIntan(fullfile(dataFold,'Ephys'),animalId,mergeId,probe,'BRN')

%Load merge info
load(fullfile(dataFold,'Ephys',animalId,mergeName,[mergeName '_mergeInfo.mat']))
nFiles = length(mergeInfo.files);
goodUnits = [];
for f = 1:nFiles
    exptName{f,1} = [animalId '_' mergeInfo.files{f}];

    sumStats{f} = anaOri(animalId,exptName{f}(8:10),exptName{f}(12:14),probe,anaMode,dataFold,0,0,f);
    goodUnits = [goodUnits screenUnits(sumStats{f},anaMode)];
end
% %remove units that do not pass inclusion criteria in either control or cool
% goodIdx = goodUnits(:,1)|goodUnits(:,2);
% for f = 1:nFiles
%     sumStats{f} = sumStats{f}(goodIdx,:);
% end

%% Plot 
if plt == 1

    clrs = {'k','c','k'};
    linStyls = {'-','-','--'};
    plr = 0;
    alignTC = 0;

    nU = height(sumStats{1});
    for u = 1:nU
        figure; hold on
        nr = 2; nc = 2;
        trialCount = 0;
        for f = 1:nFiles
    
            dat = sumStats{f};
            clr = clrs{f};
            linStyl = linStyls{f};
    
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
            subplot(nr,nc,1);hold on; box on
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
            subplot(nr,nc,3);hold on; box on
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
    
    
    %RASTER (condition)
%             subplot(nr,nc,4);hold on
%             for t = 1:nTrials
%                 xT = x(y==t & idx);
%                 if isempty(xT)
%                     continue
%                 end
%                 yT = unique(dat.fr(u).trialCond(dat.fr(u).trialNum==t));
%                 plot(xT,yT,'.','Color',clr)
%             end
%             xlim([-predelay stimTime+postdelay])
%             if ~isempty(dat.rBlank{u})
%                 ylim([0 nConds+1])
%             else
%                 ylim([0 nConds])
%             end
%             yticIdx = 1:2:nConds;
%             yticks(yticIdx)
%             yticklabels(num2str(dat.cndKey{u}(yticIdx,:)))
%             patch([0 1 1 0],[0 0 nConds+1 nConds+1],'k','EdgeColor','none','FaceAlpha',0.2)
%             xlabel('time (sec)')
%             ylabel('condition')
    
    %TUNING CURVE
%             if f==3
%                 continue
%             end
            x = dat.condition{u}(strcmp(dat.paramKey{u},'ori'),:);
            y = dat.response{u};
            meanY = mean(y,'omitnan');
            x_pref = dat.oriPref(u);
            y_pref = dat.rPref(u);
            sem = std(y,'omitnan')/sqrt(size(y,1));
            if plr == 1
                subplot(nr,nc,2,polaraxes);hold on
        %         polarplot(deg2rad(x),y,'.','Color',clr)
                polarplot(deg2rad([x x(1)]),mean([y y(:,1)],'omitnan') ,'-o','Color',clr)
                polarplot(repmat(deg2rad(x),2,1),mean(y,'omitnan')+([1;-1]*sem),'Color',clr)
                polarplot(deg2rad(dat.oriPref(u)),dat.rPref(u),'r*')
                polarplot(deg2rad([0 dat.meanVec{u}.angDir]),[0 dat.meanVec{u}.magDir],'k','LineWidth',2)
                polarplot(deg2rad([0 dat.meanVec{u}.angDir]),[0 dat.meanVec{u}.magDir*dat.ldr(u)],'g','LineWidth',2)
            else
                subplot(nr,nc,2);hold on;box on
                if alignTC == 1
                    [x,meanY,i] = alignDirTuning(x,meanY);
                    sem = sem(i);
                    y = y(:,i);
                else
                    x = [x 360];
                    y = [y y(:,1)];
                    meanY = [meanY meanY(1)];
                    sem = [sem sem(1)];
                end
%                 plot(x,y,'.','Color',clr)
                plot(x,meanY,'-square','Color',clr,'MarkerSize',7,'MarkerFaceColor',clr,'LineStyle',linStyl,'LineWidth',2)
                plot(repmat(x,2,1),meanY+([1;-1]*sem),'Color',clr,'LineWidth',2)
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
    
            trialCount = trialCount+nTrials;
        end
        ttl = [dat.exptName{u}(1:5) ' ' mergeId ' ' dat.area{u} ' ' dat.uInfo{u} ' ' num2str(dat.uID(u))];
        sgtitle(ttl)
    end

end




end