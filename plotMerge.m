%plot merged
clear all
close all

if ispc
    dataFold = 'D:\data'; 
%     dataFold = 'C:\Users\brand\Documents\data';
%     dataFold = 'F:\Brandon\data';
elseif ismac
%     dataFold = '/Volumes/Lab drive/Brandon/data';
    dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end

animal = 'febg7';
units = {'000','000','000','000','000','000'};
expts = {'000','003','004','008','009','012'};
grp = [1 1 2 2 3 3];
mergeID = '000000000003000004000008000009000012';
% units = {'001','001','001','001','001','001'};
% expts = {'000','004','005','008','009','012'};
% grp = [1 1 2 2 3 3];
% mergeID = '001000001004001005001008001009001012';




% animal = 'febg8';
% units = {'001','001','001','001','001','001'};
% expts = {'000','003','004','007','008','011'};
% grp = [1 1 2 2 3 3];
% mergeID = '001000001003001004001007001008001011';




% animal = 'febg9';
% % units = {'000','000','000','000','000','000'};
% % expts = {'002','006','007','010','011','014'};
% % grp =   [    1,    1,    2,    2,    3,    3];
% % mergeID = '000002000006000007000010000011000014';
% units = {'000','000','000'};
% expts = {'006','007','011'};
% grp =   [    1,    2,    3];
% mergeID = '000006000007000011';



% animal = 'febl0';
% % units = {'001','001','001'};
% % expts = {'006','010','016'};
% % grp = [1 2 3];
% % mergeID = '001006001010001016';
% 
% % units = {'001','001'};
% % expts = {'018','019'};
% % grp = [1 2];
% % mergeID = '001018001019';
% 
% units = {'000','000','000'};
% expts = {'010','012','015'};
% grp = [1 2 3];
% mergeID = '000010000012000015';



% animal = 'febh5';
% % units = {'000','000','000'};
% % expts = {'003','006','009'};
% % grp =   [1 2 3];
% % mergeID = '000003000006000009';
% units = {'001','001','001'};
% expts = {'004','008','015'};
% grp =   [1 2 3];
% mergeID = '001004001008001015';






% animal = 'febf6';
% units = {'001','001','001'};
% expts = {'000','003','007'};
% grp =   [1 2 3];
% mergeID = '001000001003001007';


% animal = 'febj8';
% units = {'003','003','003','003','003','003','003','003','003','003','003','003','003'};
% expts = {'002','003','004','005','006','008','009','010','016','017','018','019','020'};
% grp =   [    1,    1,    1,    1,    1,    2,    2,    2,    3,    3,    3,    3,    3];
% mergeID = '003002003003003004003005003006003008003009003010003016003017003018003019003020';
% % units = {'003','003','003','003','003','003','003','003'};
% % expts = {'002','003','004','005','006','008','009','010'};
% % grp =   [    1,    1,    1,    1,    1,    2,    2,    2];
% % mergeID = '003002003003003004003005003006003008003009003010';





% animal = 'febg2';
% units = {'001','001','001','001','001'};
% expts = {'013','016','017','020','021'};
% grp =   [1 2 3 2 3];
% mergeID = '001013001016001017001020001021';



% animal = 'febg3';
% units = {'001','001','001'};
% expts = {'011','013','017'};
% grp =   [1 2 3];
% mergeID = '001011001013001017';




clr = {'k','c','m'};
mergeName = [animal '_uMMM_' mergeID];
physDir = fullfile(dataFold,'Ephys');

binWidth=0.010; %sec
% 1sec baseline and 1sec stim period
startBin=ceil(-1/binWidth)*binWidth; %need multiple of binWidth to make 0 an edge
stopBin=floor(1/binWidth)*binWidth;
binVec=[startBin:binWidth:stopBin];
anaMode = 'SU';


load(fullfile(physDir,animal,mergeName,[mergeName '_id.mat']))
for p = 1:length(id.probes)
    if strcmp(id.probes(p).area,'PSS')
        probe = p;
    end
end
if strcmp(anaMode,'SU')
    load(fullfile(physDir,animal,mergeName,[mergeName '_p' num2str(probe) '_spkSort.mat']))
    unitIDs = unique(spkSort.unitid(spkSort.unitid~=0));
    clear spkSort
elseif strcmp(anaMode,'MU')
    unitIDs = 1:id.probes(probe).nChannels;
end


countU = 0;
x_lb = 0;
x_ub = 0;
for f = 1:length(expts)

    exptName{f} = [animal '_u' units{f} '_' expts{f}];
    load(fullfile(physDir,animal,exptName{f},[exptName{f} '_id.mat']),'id')
    for p = 1:length(id.probes)
        if strcmp(id.probes(p).area,'PSS')
            probe = p;
        end
    end
    load(fullfile(physDir,animal,exptName{f},[exptName{f} '_trialInfo.mat']),'trialInfo')
    load(fullfile(physDir,animal,exptName{f},[exptName{f} '.analyzer']),'-mat')
    if strcmp(anaMode,'SU')
        load(fullfile(physDir,animal,mergeName,[exptName{f} '_p' num2str(probe) '_' num2str(f) '_spkSort.mat']),'spkSort')
        spkStrct = spkSort;
        clear spkSort
    elseif strcmp(anaMode,'MU')
        load(fullfile(physDir,animal,exptName{f},[exptName{f} '_p' num2str(probe) '_MUspkMerge.mat']),'MUspkMerge')
        spkStrct = MUspkMerge;
        clear MUspkMerge
    end
    
    sf = id.sampleFreq;
    % kernel = ones(1,0.1*sf)*(1/(0.1*sf));
    kernel = normpdf(-3:6/2000:3);
    colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
    nTrials = length(trialInfo.triallist);
    nConds = length(unique(trialInfo.triallist));
    nDom = length(trialInfo.dom);
    nReps = nTrials/nConds;
    if strcmp(anaMode,'SU')
        uIDs = unique(spkStrct.unitid); uIDs = uIDs(uIDs~=0);
    elseif strcmp(anaMode,'MU')
        uIDs = unique(spkStrct.detCh);
    end

    predelay = getparam('predelay',Analyzer); x_lb = max(predelay,x_lb);
    stimTime = getparam('stim_time',Analyzer);
    postdelay = getparam('postdelay',Analyzer); x_ub = max(stimTime+postdelay,x_ub);
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

    for u = unitIDs

        stimCent{f,u} = [];
        countU = countU+1;
        fileID{countU,1} = [f grp(f) u];
        exptID{countU,1} = exptName{f};
        uID(countU) = u;

        if ~ismember(u,uIDs)

            info{countU,1} = 'not detected';
            raster{countU,1} = stimCent{f,u};
            fr{countU,1} = struct(  'trial',sortTrialInd,...
                                    'condition',sortTrialCond,...
                                    'stim',zeros(nReps,nConds),...
                                    'base',zeros(nReps,nConds),...
                                    'bcfr',zeros(nReps,nConds)   );
            x{countU,1} = trialInfo.domval(sortTrialCond(:,~blank));
            y{countU,1} = zeros(nReps,sum(~blank));
            rBlank{countU,1} = zeros(nReps,sum(blank));
            rPref(countU) = 0;
            cPref(countU) = nan;
            cNull(countU) = nan;
            rNull(countU) = 0;
            dsi(countU) = nan;
            dcv(countU) = nan;
            goodUnit(countU) = false;
            latCh(countU)= nan;
            latCh2(countU)= nan;

        else

            if strcmp(anaMode,'SU')
                info{countU,1} = spkStrct.unitinfo{u};
                spkTimes = spkStrct.spktimes(spkStrct.unitid == u);
            elseif strcmp(anaMode,'MU')
                info{countU,1} = 'MU';
                spkTimes = spkStrct.spktimes(spkStrct.detCh == u);
            end
    
            for c = 1:nConds
                trials = find(trialInfo.triallist == c);
                for r = 1:nReps
                    t = trials(r);
                    tvPre       = stimStart(t)-(predelay*sf):stimStart(t)-1;
                    tvStim      = stimStart(t):stimStart(t)+(stimTime*sf)-1;
                    tvPost      = stimStart(t)+(stimTime*sf):stimStart(t)+((stimTime+postdelay)*sf)-1;
                    tvTrial     = [tvPre tvStim tvPost];
                    spkTrain{f}(:,t,u) = ismember(tvTrial,spkTimes);
    
                    stimCent{f,u} = [stimCent{f,u} vertcat( (tvTrial(spkTrain{f}(:,t,u))-stimStart(t))/sf , ...
                                                                    repmat(t,1,length(find(spkTrain{f}(:,t,u)))) ) ];
    
                    baseCount = length(find(ismember(stimStart(t)-sf:stimStart(t),spkTimes)));
                    stimCount = length(find(ismember(stimStart(t):stimStart(t)+sf,spkTimes)));
                    baseFR{f}(r,c,u) = baseCount;
                    stimFR{f}(r,c,u) = stimCount;
                    bcfr{f}(r,c,u) = stimCount-baseCount;
    
                end
    
            end
            raster{countU,1} = stimCent{f,u};
            fr{countU,1} = struct(  'trial',sortTrialInd,...
                                    'condition',sortTrialCond,...
                                    'stim',stimFR{f}(:,:,u),...
                                    'base',baseFR{f}(:,:,u),...
                                    'bcfr',bcfr{f}(:,:,u)   );
            x{countU,1} = trialInfo.domval(sortTrialCond(:,~blank));
            y{countU,1} = bcfr{f}(:,~blank,u);
            rBlank{countU,1} = stimFR{f}(:,blank,u);
            rPref(countU) = max(mean(y{countU,1}));
            cP = x{countU,1}(1,mean(y{countU,1}) == rPref(countU));
            if length(cP)>1 % if there is more than one pk with rPref
                rIn = mean(y{countU,1},'omitnan');
                pks = rIn == rPref(countU);
                rConv = rIn([end 1:end 1]);
                rConv = conv(rConv,ones(1,3)*(1/3),'same');
                rConv = rConv(1+1:end-1);
                rConv(~pks) = 0;
                cPref(countU) = x{countU,1}(1,find(rConv==max(rConv),1,'first'));
                clear rIn pks rConv 
            elseif length(cP)==1
                cPref(countU) = cP;
            end
            cNull(countU) = mod(cPref(countU)+180,360);
            rNull(countU) = mean(y{countU,1}( : , mean(x{countU,1})==cNull(countU) ));
            dsi(countU) = abs(rPref(countU)-rNull(countU)) / rPref(countU);
            mv = meanvec(mean(x{countU,1}),mean(y{countU,1}));
            dcv(countU) = mv.cv;
    %         [g] = dirGauss(mean(y),mean(x),0);
    %         xMdl = linspace(0,359,360);
            isAct = rPref(countU)>=2;
            bfr = fr{countU,1}.base(:,~blank); bfr = bfr(:);
            sfr = fr{countU,1}.stim(:,~blank); sfr = sfr(:);
            pVis = ranksum(bfr,sfr);
            if isnan(pVis)
                pVis = 1;
            end
            isVis = pVis<0.01;
    
            if strcmp(anaMode,'SU')
                isSU = strcmp(info{countU,1},'SU');
                goodUnit(countU) = isAct & isVis & isSU;
            else
                goodUnit(countU) = isAct & isVis;
            end
    
            prefTrials = find(trialInfo.triallist== find(trialInfo.domval == cPref(countU)) );
            for rep = 1:length(prefTrials)
                t = prefTrials(rep);
                spkTs{rep} = raster{countU,1}(1,raster{countU,1}(2,:)==t);
                N(rep,:) = histcounts(spkTs{rep},binVec);
            end
            avgN = mean(N,1)/binWidth;
            avgBase=mean(avgN(binVec<0));
            stdBase=std(avgN(binVec<0));
            cumN=cumsum(avgN-avgBase); %cumsum(1): value 1 in input
            diffN=diff(cumN);
            %criterion 1: cumsum over threshold and 2 increasing bins
            idx=find(cumN(1:end-2)>2*stdBase & diffN(1:end-1)>0 & diffN(2:end)>0 ...
                & binVec(1:end-3)>0,1);
            if ~isempty(idx)
                latCh(countU)=binVec(idx);
            else
                latCh(countU)=NaN;
            end
    
            %criterion 2: cumsum over threshold and 3 increasing bins
            idx=find(cumN(1:end-3)>2*stdBase & ...
                diffN(1:end-2)>0 & diffN(2:end-1)>0 & diffN(3:end)>0 ...
                & binVec(1:end-4)>0,1);
            if ~isempty(idx)
                latCh2(countU)=binVec(idx);
            else
                latCh2(countU)=NaN;
            end

        end

    end

end

varNames = {'fileID','exptID','uID','uInfo','goodUnit','raster','lat1','lat2','fr','tuningX','tuningY','rBlank','rPref','cPref','rNull','DSI','DCV'};
uDat = table(fileID,exptID,uID',info,goodUnit',raster,latCh',latCh2',fr,x,y,rBlank,rPref',cPref',rNull',dsi',dcv','VariableNames',varNames);

clear x y uIDs uID exptID h

uIDs = unique(uDat.uID);
for g = unique(grp)
    gInd = ismember(uDat.exptID,exptName(grp == g));
    GU(:,g) = ismember(uIDs,uDat.uID(gInd & uDat.goodUnit));
end
keepUnits = uIDs(GU(:,1)|GU(:,2));
uDat = uDat(ismember(uDat.uID,keepUnits),:);
uIDs = keepUnits;
for u = 1:length(uIDs)
    figure;hold on
    countTrial = 0;
    ttl = [];
    swtch = 1;
    for g = unique(grp)

        curDat = uDat(ismember(uDat.exptID,exptName(grp==g)') & uDat.uID==uIDs(u),:);
        if height(curDat) == 0
            continue
        end

        subplot(2,2,4);hold on
        for e = 1:height(curDat)
            plot(mean(curDat.tuningX{e}),mean(curDat.tuningY{e}),'--','Color',clr{g})
            plot(repmat(mean(curDat.tuningX{e}),2,1),mean(curDat.tuningY{e})+([-1;1]*sem(curDat.tuningY{e})),'Color',clr{g})
        end
        x = vertcat(curDat.tuningX{:});
        y = vertcat(curDat.tuningY{:});
        plot(x(:),y(:),'.','Color',clr{g})
        plot(curDat.cPref,curDat.rPref,'*','Color',clr{g})

        subplot(2,2,2);hold on
        plot(mean(x),mean(y),'LineWidth',2,'Color',clr{g})
        plot(repmat(mean(x),2,1),mean(y)+([-1;1]*sem(y)),'LineWidth',2,'Color',clr{g})
%         ttl = [ttl 'DSI(' clr{g} ') = ' num2str(mean(curDat.DSI)) ';'];
%         title(ttl)

        subplot(2,2,3);hold on
        for f = 1:height(curDat)
            nTrial = max(max(curDat.fr{f}.trial));
            if isempty(curDat.raster{f})
                countTrial = countTrial+nTrial;
                continue
            end
            
            x = curDat.raster{f}(1,:);
            y = curDat.raster{f}(2,:)+countTrial;
            countTrial = countTrial+nTrial;
            plot(x,y,'.','Color',clr{g})
        end
        xlim([-1*x_lb x_ub])
        ylim([0 countTrial])

        subplot(2,2,1);hold on
        h = [curDat.raster{:}];
        if ~isempty(h)
            h = h(1,:);
        end
        binW = 0.1;
        hist = histogram(h,'BinWidth',binW,'FaceColor',clr{g},'FaceAlpha',0.3,'EdgeColor','none');
        hist.BinCounts = (hist.BinCounts/nTrial)/binW;
        xlim([-1*x_lb x_ub])

        if swtch == 1
            uTypeTemp = unique(vertcat(curDat.uInfo(curDat.uID==uIDs(u)))); uType{uIDs(u)} = uTypeTemp{:};
            sgtitle([animal ': unit#' num2str(uIDs(u)) ' (' uType{uIDs(u)} ')'])
            swtch = 0;
        end

    end

    saveas(gcf,fullfile(physDir,animal,mergeName,[anaMode num2str(uIDs(u)) 'plot']),'fig')
end


figure;hold on
for g = unique(grp)
    curDat = uDat( ismember(uDat.exptID,exptName(grp==g)') ,:);
    for u = unitIDs
        if ismember(u,curDat.uID)
            LAT(g,u) = mean(curDat(curDat.uID == u,:).lat1,'omitnan');
            RP(g,u) = mean(curDat(curDat.uID == u,:).rPref);
            RN(g,u) = mean(curDat(curDat.uID == u,:).rNull);
            RB(g,u) = mean(mean([curDat(curDat.uID == u,:).rBlank{:}],1),2);
            DSI(g,u) = mean(curDat(curDat.uID == u,:).DSI,'omitnan');
            DCV(g,u) = mean(curDat(curDat.uID == u,:).DCV,'omitnan');
        else
            LAT(g,u) = nan;
            RP(g,u) = nan;
            RN(g,u) = nan;
            RB(g,u) = nan;
            DSI(g,u) = nan;
            DCV(g,u) = nan;
        end
    end

    subplot(2,2,1);hold on
    if ~( sum(isnan(LAT(g,:))) == length(LAT(g,:)) )
        cdf = cdfplot(LAT(g,:));
        cdf.Color = clr{g};
    end
    title('')
    ylabel('percentile')
    xlabel('latency')

    subplot(2,2,2);hold on
    cdf = cdfplot(RP(g,:));
    cdf.Color = clr{g};
    title('')
    ylabel('percentile')
    xlabel('rPref')

    subplot(2,2,3);hold on
    cdf = cdfplot(RB(g,:));
    cdf.Color = clr{g};
    title('')
    ylabel('percentile')
    xlabel('rBlank')

    subplot(2,2,4);hold on
    cdf = cdfplot(1-DCV(g,:));
    cdf.Color = clr{g};
    title('')
    ylabel('percentile')
    xlabel('1-DCV')
    xlim([0 1])

end

figure; hold on
subplot(2,2,1);hold on
x = LAT(1,:);
y = LAT(2,:);
plot(x,y,'k.')
ub = max([max(x) max(y)]);
plot([0 ub],[0 ub],'k--')
xlabel('latency before cooling V1')
ylabel('latency during cooling V1')

subplot(2,2,2);hold on
x = RP(1,:);
y = RP(2,:);
plot(x,y,'k.')
ub = max([max(x) max(y)]);
plot([0 ub],[0 ub],'k--')
xlabel('rPref before cooling V1')
ylabel('rPref during cooling V1')

subplot(2,2,3);hold on
x = RB(1,:);
y = RB(2,:);
plot(x,y,'k.')
ub = max([max(x) max(y)]);
plot([0 ub],[0 ub],'k--')
xlabel('rBlank before cooling V1')
ylabel('rBlank during cooling V1')

subplot(2,2,4);hold on
x = 1-DCV(1,:);
y = 1-DCV(2,:);
plot(x,y,'k.')
ub = max([max(x) max(y)]);
plot([0 ub],[0 ub],'k--')
xlabel('1-DCV before cooling V1')
ylabel('1-DCV during cooling V1')

sgtitle([animal ' summary plot'])
saveas(gcf,fullfile(physDir,animal,mergeName,[anaMode 'summaryPlot']),'fig')

save(fullfile(physDir,animal,mergeName,[anaMode 'dat.mat']),'LAT','RP','RB','DSI','DCV','GU','uDat')
