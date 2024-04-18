%plot merged
clear all
close all

if ispc
    dataFold = 'D:\data'; 
%     dataFold = 'F:\Brandon\data';
elseif ismac
%     dataFold = '/Volumes/Lab drive/Brandon/data';
    dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end
% animal = 'febj8';
% units = {'003','003','003','003','003','003','003','003','003','003','003','003','003'};
% expts = {'002','003','004','005','006','008','009','010','016','017','018','019','020'};
% grp =   [    1,    1,    1,    1,    1,    2,    2,    2,    3,    3,    3,    3,    3];
animal = 'febg9';
units = {'000','000','000','000','000','000'};
expts = {'002','006','007','010','011','014'};
grp =   [    1,    1,    2,    2,    3,    3];
clr = {'k','c','m'};
probe = 1;
mergeID = [];
for f = 1:length(expts)
    mergeID = [mergeID units{f}];
    mergeID = [mergeID expts{f}];
end
mergeName = [animal '_uMMM_' mergeID];
physDir = fullfile(dataFold,'Ephys');

countU = 0;
for f = 1:length(expts)

    exptName{f} = [animal '_u' units{f} '_' expts{f}];
    load(fullfile(physDir,animal,exptName{f},[exptName{f} '_id.mat']),'id')
    load(fullfile(physDir,animal,exptName{f},[exptName{f} '_trialInfo.mat']),'trialInfo')
    load(fullfile(physDir,animal,exptName{f},[exptName{f} '.analyzer']),'-mat')
    load(fullfile(physDir,animal,mergeName,[exptName{f} '_p' num2str(probe) '_' num2str(f) '_spkSort.mat']),'spkSort')
    
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
        stimCent{f,u} = [];
        countU = countU+1;
        exptID{countU,1} = exptName{f};
        uID(countU) = uIDs(u);
        info{countU,1} = spkSort.unitinfo{uIDs(u)};
        spkTimes = spkSort.spktimes(spkSort.unitid == uIDs(u));
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
        rBlank{countU,1} = bcfr{f}(:,blank,u);
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
        rNull(countU) = mean(y{countU,1}(:,mean(x{countU,1})==cNull(countU)));
        dsi(countU) = (rPref(countU)-rNull(countU))/rPref(countU);
%         [g] = dirGauss(mean(y),mean(x),0);
%         xMdl = linspace(0,359,360);

    end

end

varNames = {'exptID','uID','uInfo','raster','fr','tuningX','tuningY','rBlank','rPref','cPref','DSI'};
uDat = table(exptID,uID',info,raster,fr,x,y,rBlank,rPref',cPref',dsi','VariableNames',varNames);

clear x y uIDs uID exptID

uIDs = unique(uDat.uID);
for u = 1:length(uIDs)
    figure;hold on
    countTrial = 0;
    for g = unique(grp)

        curDat = uDat(ismember(uDat.exptID,exptName(grp==g)') & uDat.uID==uIDs(u),:);
        
        subplot(2,2,4);hold on
        for e = 1:height(curDat)
            plot(mean(curDat.tuningX{e}),mean(curDat.tuningY{e}),'--','Color',clr{g})
            plot(repmat(mean(curDat.tuningX{e}),2,1),mean(curDat.tuningY{e})+([-1;1]*sem(curDat.tuningY{e})),'Color',clr{g})
        end
        x = vertcat(curDat.tuningX{:});
        y = vertcat(curDat.tuningY{:});
        plot(x(:),y(:),'.','Color',clr{g})

        subplot(2,2,2);hold on
        plot(mean(x),mean(y),'LineWidth',2,'Color',clr{g})
        plot(repmat(mean(x),2,1),mean(y)+([-1;1]*sem(y)),'LineWidth',2,'Color',clr{g})

        subplot(2,2,3);hold on
        for f = 1:height(curDat)
            nTrial = size(curDat.tuningX{f},1)*size(curDat.tuningX{f},2)+length(curDat.rBlank{f});
            x = curDat.raster{f}(1,:);
            y = curDat.raster{f}(2,:)+countTrial;
            countTrial = countTrial+nTrial;
            plot(x,y,'.','Color',clr{g})
        end

        subplot(2,2,1);hold on
        h = [curDat.raster{:}];h = h(1,:);
        histogram(h,'BinWidth',.1,'FaceColor',clr{g},'FaceAlpha',0.3,'EdgeColor','none')

    end
    sgtitle(['unit#' num2str(uIDs(u))])

    saveas(gcf,fullfile(physDir,animal,mergeName,['u' num2str(u) 'plot']),'fig')
end

