
clear all
close all

%Settings
animalId = 'febj7';
mergeId = '000015000017000018';
probeId = 1;
mergeName = [animalId '_uMMM_' mergeId];
dataFold = 'F:\Brandon\data';
anaMode = 'SU';

clrs = {'k','c','m'};
plr = 0;
alignTC = 0;

%Load merge info
load(fullfile(dataFold,'Ephys',animalId,mergeName,[mergeName '_mergeInfo.mat']))
nFiles = length(mergeInfo.files);
goodUnits = [];
for f = 1:nFiles
    exptName{f,1} = [animalId '_' mergeInfo.files{f}];

    sumStats{f} = anaOri(animalId,exptName{f}(8:10),exptName{f}(12:14),probeId,anaMode,dataFold,0,0,f);
    goodUnits = [goodUnits sumStats{f}.goodUnit];
end

nU = height(sumStats{1});
for u = find(goodUnits(:,1)|goodUnits(:,2))'
    figure; hold on
    nr = 2; nc = 2;
    trialCount = 0;
    for f = 1:nFiles

        dat = sumStats{f};
        clr = clrs{f};

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

        subplot(nr,nc,1);hold on
        binSize = 0.1;
        bins = -predelay:binSize:stimTime+postdelay;
        h = histogram(x,'BinEdges',bins);
        h.EdgeColor = 'none';
        h.FaceColor = clr;
        h.BinCounts = h.BinCounts/(nTrials*binSize);
        if ~isnan(dat.latency(u))
            xline(dat.latency(u),'--','LineWidth',2)
        end
        xlim([bins(1) bins(end)])


        subplot(nr,nc,3);hold on
        idx = x>-predelay & x<(stimTime+postdelay);
        plot(x(idx),y(idx),'.','Color',clr)
        trialExclude = dat.fr(u).trialNum(  isnan(dat.fr(u).bc)  );
        for t = trialExclude'
            patch([-predelay stimTime+postdelay stimTime+postdelay -predelay],[t-0.5 t-0.5 t+0.5 t+0.5],'r','EdgeColor','none','FaceAlpha',0.2)
        end
        if ~isnan(dat.latency(u))
            xline(dat.latency(u),'--','LineWidth',2)
        end
        xlim([-predelay stimTime+postdelay])
        ylim([0 nTrials+trialCount])
        patch([0 1 1 0],[trialCount trialCount nTrials nTrials],'k','EdgeColor','none','FaceAlpha',0.2)
        xlabel('time (sec)')
        ylabel('trial #')


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
            subplot(nr,nc,2);hold on
            if alignTC == 1
                [x,meanY,i] = alignDirTuning(x,meanY);
                sem = sem(i);
                y = y(:,i);
            end
            plot(x,y,'.','Color',clr)
            plot(x,meanY,'-o','Color',clr)
            plot(repmat(x,2,1),meanY+([1;-1]*sem),'Color',clr)
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

        trialCount = trialCount+nTrials;
    end
end