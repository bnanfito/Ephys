%plot merged
clear all
close all

dataFold = 'C:\Users\brand\Documents\data';
animal = 'febg9';
units = {'000','000'};
expts = {'006','011'};
probe = 1;
mergeID = [];
for f = 1:length(expts)
    mergeID = [mergeID units{f}];
    mergeID = [mergeID expts{f}];
end
mergeName = [animal '_uMMM_' mergeID];

physDir = fullfile(dataFold,'Ephys');
clr = {'k','m'};
for f = 1:length(expts)

    exptName = [animal '_u' units{f} '_' expts{f}];
    load(fullfile(physDir,animal,exptName,[exptName '_id.mat']),'id')
    load(fullfile(physDir,animal,exptName,[exptName '_trialInfo.mat']),'trialInfo')
    load(fullfile(physDir,animal,exptName,[exptName '.analyzer']),'-mat')
    load(fullfile(physDir,animal,mergeName,[exptName '_p' num2str(probe) '_' num2str(f) '_spkSort.mat']),'spkSort')
    
    area = id.probes(probe).area;
    sf = id.sampleFreq;
    % kernel = ones(1,0.1*sf)*(1/(0.1*sf));
    kernel = normpdf(-3:6/2000:3);
    colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
    nTrials = length(trialInfo.triallist);
    nConds = length(unique(trialInfo.triallist));
    nDom = length(trialInfo.dom);
    nReps = nTrials/nConds;
    nU = length(spkSort.unitinfo);
    uIDs = unique(spkSort.unitid); uIDs = uIDs(uIDs~=0);

    predelay = getparam('predelay',Analyzer);
    stimTime = getparam('stim_time',Analyzer);
    postdelay = getparam('postdelay',Analyzer);
    if stimTime < 1
        st = stimTime;
    else
        st = 1;
    end
    trialStart = trialInfo.eventTimes([1;find(diff(trialInfo.eventId) == 1)+1]);
    stimStart = trialInfo.eventTimes(find(diff(trialInfo.eventId) == 2)+1);
    stimEnd = trialInfo.eventTimes(find(diff(trialInfo.eventId) == -2)+1);
    trialEnd = trialInfo.eventTimes(find(diff(trialInfo.eventId) == -1)+1);
    [sortTrialCond,sortTrialInd] = sort(trialInfo.triallist);
    sortTrialCond = reshape(sortTrialCond,nReps,nConds);
    sortTrialInd = reshape(sortTrialInd,nReps,nConds);
    trialL = predelay+stimTime+postdelay;
    nT = ((1:(trialL*sf))-(predelay*sf))/sf;
    if isempty(trialInfo.blankId)
        blank = zeros(1,nConds)==1;
    else
        blank = (1:nConds)==trialInfo.blankId;
    end

    for u = 1:length(uIDs)
        uID = uIDs(u);
        figure(u);hold on
        spkTimes = spkSort.spktimes(spkSort.unitid == uID);
        for c = 1:nConds
            trials = find(trialInfo.triallist == c);
            for r = 1:nReps
                t = trials(r);
                tvPre       = stimStart(t)-(predelay*sf):stimStart(t)-1;
                tvStim      = stimStart(t):stimStart(t)+(stimTime*sf)-1;
                tvPost      = stimStart(t)+(stimTime*sf):stimStart(t)+((stimTime+postdelay)*sf)-1;
                tvTrial     = [tvPre tvStim tvPost];
                spkTrain{f}(:,u) = ismember(tvTrial,spkTimes);
                baseCount = length(find(ismember(stimStart(t)-sf:stimStart(t),spkTimes)));
                stimCount = length(find(ismember(stimStart(t):stimStart(t)+sf,spkTimes)));
                baseFR{f}(r,c,u) = baseCount;
                stimFR{f}(r,c,u) = stimCount;
                bcfr{f}(r,c,u) = stimCount-baseCount;

            end

        end
        x = trialInfo.domval(sortTrialCond(:,~blank));
        y = bcfr{f}(:,~blank,u);
        rBlank = bcfr{f}(:,blank,u);
        rPref = max(mean(y));
        cPref = x(1,mean(y) == rPref);
        if length(cPref)>1 % if there is more than one pk with rPref
            rIn = mean(y,'omitnan');
            pks = rIn == rPref;
            rConv = rIn([end 1:end 1]);
            rConv = conv(rConv,ones(1,3)*(1/3),'same');
            rConv = rConv(1+1:end-1);
            rConv(~pks) = 0;
            cPref = x(1,find(rConv==max(rConv),1,'first'));
            clear rIn pks rConv 
        end
        cNull = mod(cPref+180,360);
        rNull = mean(y(:,mean(x)==cNull));
        dsi = (rPref-rNull)/rPref;
        [g] = dirGauss(mean(y),mean(x),0);
        xMdl = linspace(0,359,360);


%         plot(x,y,[clr{f} '.'],'MarkerSize',10)
        patch([x(1,1) x(1,end) x(1,end) x(1,1)],[min(rBlank) min(rBlank) max(rBlank) max(rBlank)],...
            'r','EdgeColor','none','FaceAlpha',0.2)
        plot(xMdl,g.auss(xMdl),'-','LineWidth',2,'Color',clr{f})
        plot(mean(x),mean(y),'--o','LineWidth',1,'Color',clr{f})
        plot(repmat(mean(x),2,1) , mean(y)+([-1;1]*(std(y)/sqrt(size(y,1)))) ,'-','LineWidth',2,'Color',clr{f})
        plot(cPref,rPref,'*','MarkerSize',10,'LineWidth',2,'Color',clr{f})
        xlabel('dir. (deg.)')
        ylabel('BCFR (Hz)')
        yline(0,'k--')
        title(['unit #' num2str(u)])

        uDat{f}(u).uID = uID;
        uDat{f}(u).fr = struct(  'trial',sortTrialInd,...
                                'condition',sortTrialCond,...
                                'stim',stimFR{f}(:,:,u),...
                                'base',baseFR{f}(:,:,u),...
                                'bcfr',bcfr{f}(:,:,u)   );
        uDat{f}(u).rPref = rPref;
        uDat{f}(u).cPref = cPref;
        uDat{f}(u).tCurveX = mean(x);
        uDat{f}(u).tCurveY = y;
        uDat{f}(u).rBlank = rBlank;
        uDat{f}(u).dsi = dsi;

    end


end

figure;

subplot(2,2,1);hold on
plot([uDat{1}.rPref],[uDat{2}.rPref],'k.')
plot([0 20],[0 20],'k--')
xlabel('rPref Before (Hz)')
ylabel('rPref During (Hz)')

subplot(2,2,2);hold on
plot([uDat{1}.cPref],[uDat{2}.cPref],'k.')
plot([0 360],[0 360],'k--')
xlabel('cPref Before (deg)')
ylabel('cPref During (deg)')

subplot(2,2,3);hold on
plot([uDat{1}.dsi],[uDat{2}.dsi],'k.')
plot([0 1],[0 1],'k--')
xlabel('dsi Before')
ylabel('dsi During')





