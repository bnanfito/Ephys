% anaOri
function [sumStats] = anaOri(animal,unit,expt,probe,anaMode,dataFold,plt,svePlt)

%% Initialize

% clear
% close all

% if ispc
% %     dataFold = 'D:\data'; 
%     dataFold = 'C:\Users\brand\Documents\data';
% %     dataFold = 'F:\Brandon\data';
% elseif ismac
%     dataFold = '/Volumes/Lab drive/Brandon/data';
% %     dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
% end
physDir = fullfile(dataFold,'Ephys');
figDir = fullfile(dataFold,'Figures');


%% Settings

% animal = 'febl0';
% unit = '000';
% expt = '000';
exptName = [animal '_u' unit '_' expt];
% probe = 'PSS';

% plt = 1;
plr = 0;
alignBit = 0;
% anaMode = 'MU';
% stimMode = 'bi ff';
visTest = 'ranksum';
alpha = 0.01;


%% Load Data

exptDir = fullfile(physDir,animal,exptName);
load(fullfile(exptDir, [exptName '_id.mat'] ))
load(fullfile(exptDir, [exptName '_trialInfo.mat'] ))
load(fullfile(exptDir, [exptName '.analyzer']),'-mat')

varInfo = whos('probe');
if strcmp(varInfo.class,'char')
    area = probe;
    probe = find(strcmp({id.probes(:).area}',area));
end
clear varInfo
area = id.probes(probe).area;
sf = id.sampleFreq;
predelay = getparam('predelay',Analyzer);
stimTime = getparam('stim_time',Analyzer);
postdelay = getparam('postdelay',Analyzer);
trialL = predelay+stimTime+postdelay;
nTrials = length(trialInfo.triallist);
nConds = length(unique(trialInfo.triallist));
nDom = length(trialInfo.dom);
nReps = nTrials/nConds;
if isempty(trialInfo.blankId)
    blank = zeros(1,nConds)==1;
    blankTrialList = [];
else
    blank = (1:nConds)==trialInfo.blankId;
    blankTrialList = find(trialInfo.triallist==find(blank));
end
[sortTrialInd(:,1),sortTrialInd(:,2)] = sort(trialInfo.triallist);
nT = ((1:(trialL*sf))-(predelay*sf))/sf;
kernel = normpdf(-3:6/2000:3);
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};

oriInd = strcmp(trialInfo.dom,'ori');
sizeInd = strcmp(trialInfo.dom,'y_size') | strcmp(trialInfo.dom,'y_size ') | strcmp(trialInfo.dom,'x_size');
if sum(sizeInd)==0
    nSize = 1;
    szTrialList{1} = find(trialInfo.triallist~=find(blank));
elseif sum(sizeInd)>0
    sizes = unique(trialInfo.domval(:,sizeInd));
    nSize = length(sizes);
    for sz = 1:nSize
        szTrialList{sz} = find(ismember(trialInfo.triallist,find(trialInfo.domval(:,sizeInd)==sizes(sz))));
    end
end

%% Compute 

[spks,trialExclude] = orgSpks(animal,unit,expt,probe,anaMode,dataFold);
[goodUnits,pVis] = screenUnits(spks,anaMode,blank,visTest,alpha);
nU = length(spks);

dim{1} = 'response';
for d = 1:nDom
    dim{1+d} = trialInfo.dom{d};
end

for u = 1:nU % u indexes a unit (column) in structure spks

    exptID{u} = exptName;
    uID(u) = spks(u).unitId;
    probeID(u) = probe;
    areaID{u} = area;

    for sz = 1:nSize

        if nDom==1 && strcmp(trialInfo.dom{1},'ori') && sum(sizeInd)==0
    
            cndInclude = ~blank;
    
        elseif nDom==2 && sum(sizeInd)>0
            cndInclude = trialInfo.domval( : , sizeInd ) == sizes(sz);
            if sizes(sz) >= 150
                ff = sz;
            elseif sizes(sz) <= 75
                hemi = sz;
            end
        end
    
        if ~isempty(trialInfo.blankId)
            Rblank{u} = spks(u).fr.stim(:,blank);
        end
        R{u,sz} = spks(u).fr.bc(:,cndInclude);
        C{u,sz} = trialInfo.domval(cndInclude,:)';
        paramKey{u} = trialInfo.dom;
        rMean = mean(R{u,sz},'omitnan');
        rPref = max(rMean);
        if sum(rMean==rPref)>1
            rIn = rMean;
            pks = rIn == rPref;
            rConv = rIn([end 1:end 1]);
            rConv = conv(rConv,ones(1,3)*(1/3),'same');
            rConv = rConv(1+1:end-1);
            rConv(~pks) = 0;
            cPref = C{u,sz}(oriInd,find(rConv==max(rConv),1,'first'));
            clear rIn pks rConv 
        else
            cPref = C{u,sz}(oriInd,rMean==rPref);
        end
        cNull = mod(cPref+180,360);
        rNull = rMean(C{u,sz}(oriInd,:)==cNull);
    
        dsi(u,sz) = abs(rPref-rNull)/rPref;
        mv = meanvec(C{u,sz}(oriInd,:),rMean);
        ldir(u,sz) = mv.ldir;
    
        Rpref(u,sz) = rPref;
        Cpref(u,sz) = cPref;
        Rnull(u,sz) = rNull;
        Cnull(u,sz) = cNull;

    end

end

varNames = {'exptName','probe','area','uID','goodUnit','pVis','fr','response','condition','paramKey','rPref','oriPref','rNull','oriNull','rBlank','dsi','ldr'};
sumStats = table(exptID',probeID',areaID',uID',goodUnits',pVis,vertcat(spks.fr),R,C,paramKey',Rpref,Cpref,Rnull,Cnull,Rblank',dsi,ldir,'VariableNames',varNames);

if strcmp(anaMode,'MU')
    sumStats.xPos = vertcat(spks.xPos);
    sumStats.zPos = vertcat(spks.zPos);
end


%% Plot

if plt == 1

    for u = 1:nU
        figure;
        legLbl{nSize+1} = 'blank';
        x = spks(u).stimCent(1,:);
        y = spks(u).stimCent(2,:);
        blankSpkIdx = ismember(y,blankTrialList);
        for sz = 1:nSize

            clr = colors{sz};
            if nDom==2 && sum(sizeInd)>0
                if sz==ff
                    legLbl{sz} = ['full field'];
                elseif sz==hemi
                    legLbl{sz} = ['hemi field'];
                end
            else
                legLbl{sz} = ['data'];
            end

            szSpkIdx = ismember(y,szTrialList{sz});

            subplot(2,2,1);hold on
            bins = -1:0.1:2;
            h = histogram(x(szSpkIdx),'BinEdges',bins);
            h.FaceColor = colors{sz};
            h.EdgeColor = 'none';
            xlim([bins(1) bins(end)])
    
            subplot(2,2,3);hold on
            if isempty(x) || isempty(y)
                text(0,0,'no spikes')
                ttl = [exptName ' p' num2str(probe) ' (' area ') ' anaMode '#' num2str(uID(u))];
                if ~ismember(u,find(goodUnits))
                    ttl = [ttl '(BAD UNIT)'];
                end
                sgtitle(ttl)
                continue
            end
            if sz == 1
                patch([0 1 1 0],[0 0 max(y)+1 max(y)+1],'k','EdgeColor','none','FaceAlpha',0.2)
            end
            for t = 1:nTrials
                if ismember(t,find(trialExclude))
                    patch([-predelay stimTime+postdelay stimTime+postdelay -predelay],[t-0.5 t-0.5 t+0.5 t+0.5],'r','EdgeColor','none','FaceAlpha',0.2)
                end
            end
            plot(x(szSpkIdx),y(szSpkIdx),'.','Color',colors{sz})
            xlim([-1 2])
            ylim([0 max(y)+1])

            if plr == 1
                subplot(1,2,2,polaraxes);hold on
                xT{u} = deg2rad(C{u,sz}(oriInd,:)); xT{u} = [xT{u}(1:end) xT{u}(1)];
                yT{u} = mean(R{u,sz},'omitnan'); yT{u} = [yT{u}(1:end) yT{u}(1)];
                legSubset(sz) = polarplot(xT{u},yT{u},'-o','Color',clr);
            else
                subplot(1,2,2);hold on
                xT{u} = C{u,sz}(oriInd,:);
                yT{u} = R{u,sz};
                yMean{u} = mean(yT{u},'omitnan');
                sem{u} = std(yT{u},'omitnan')/sqrt(size(yT{u},1));
                if alignBit == 1
                    [xT{u},yMean{u},i] = alignDirTuning(xT{u},yMean{u});
                    yT{u} = yT{u}(:,i);
                    sem{u} = sem{u}(i);
                end
                plot(xT{u},yT{u}','.','Color',clr)
                plot(repmat(xT{u},2,1),yMean{u}+([1;-1]*sem{u}),'Color',clr)
                legSubset(sz) = plot(xT{u},yMean{u},'-o','Color',clr);
            end

        end

        subplot(2,2,1);hold on
        bins = -1:0.1:2;
        h = histogram(x(blankSpkIdx),'BinEdges',bins);
        h.FaceColor = 'k';
        h.EdgeColor = 'none';
        xlim([bins(1) bins(end)])

        subplot(2,2,3);hold on
        plot(x(blankSpkIdx),y(blankSpkIdx),'k.')

        legend(legSubset,legLbl)

        ttl = [exptName ' p' num2str(probe) ' (' area ') ' anaMode '#' num2str(uID(u))];
        if ~ismember(u,find(goodUnits))
            ttl = [ttl '(BAD UNIT)'];
        end
        sgtitle(ttl)
        if svePlt == 1 && ismember(u,find(goodUnits))
            figFileDir = fullfile(figDir,animal,exptName);
            if ~isfolder(figFileDir)
                mkdir(figFileDir)
            end
            figFileName = fullfile(figFileDir,[exptName '_p' num2str(probe) '_' anaMode num2str(uID(u))]);
            saveas(gcf,figFileName)
        end
    end

%     figure;
%     goodIdx = sumStats.goodUnit;
%     
%     subplot(2,3,1); hold on
%     x = mean(vertcat(xT{:}),'omitnan');
%     y = mean(vertcat(yMean{:}),'omitnan');
%     plot(x,y,'k--','LineWidth',2)
%     y = mean(vertcat(yMean{goodIdx}),'omitnan');
%     plot(x,y,'k-o','LineWidth',2)
% 
%     subplot(2,3,2);hold on
%     x = sumStats.dsi;
%     cdf = cdfplot(x);
%     cdf.LineWidth = 2;
%     cdf.LineStyle = '--';
%     cdf.Color = 'k';
%     if sum(goodIdx)>0
%         cdf = cdfplot(x(goodIdx));
%         cdf.LineWidth = 2;
%         cdf.LineStyle = '-';
%         cdf.Color = 'k';
%     end
%     xlabel('DSI')
%     xlim([0 1])
%     ylabel('percentile')
%     ylim([0 1])
%     title('')
% 
%     subplot(2,3,3);hold on
%     x = sumStats.ldr;
%     cdf = cdfplot(x);
%     cdf.LineWidth = 2;
%     cdf.LineStyle = '--';
%     cdf.Color = 'k';
%     if sum(goodIdx)>0
%         cdf = cdfplot(x(goodIdx));
%         cdf.LineWidth = 2;
%         cdf.LineStyle = '-';
%         cdf.Color = 'k';
%     end
%     xlabel('Ldir')
%     xlim([0 1])
%     ylabel('percentile')
%     ylim([0 1])
%     title('')
% 
%     subplot(2,3,4); hold on
%     x = sumStats.rPref;
%     cdf = cdfplot(x);
%     cdf.LineWidth = 2;
%     cdf.LineStyle = '--';
%     cdf.Color = 'k';
%     if sum(goodIdx)>0
%         cdf = cdfplot(x(goodIdx));
%         cdf.LineWidth = 2;
%         cdf.LineStyle = '-';
%         cdf.Color = 'k';
%         legend({['all units; n=' num2str(length(x))],['good units; n=' num2str(sum(goodIdx))]})
%     else
%         legend({['no good units; n=' num2str(length(x))]})
%     end
%     xlabel('rPref')
%     ylabel('percentile')
%     ylim([0 1])
% 
%     subplot(2,3,5);hold on
%     x = sumStats.dsi;
%     y = sumStats.rPref;
%     plot(x,y,'ko','MarkerSize',5);
%     plot(x(goodIdx),y(goodIdx),'k.','MarkerSize',10)
%     xlabel('DSI')
%     xlim([0 1])
%     ylabel('rPref')
% 
%     subplot(2,3,6);hold on
%     x = sumStats.ldr;
%     y = sumStats.rPref;
%     plot(x,y,'ko','MarkerSize',5);
%     plot(x(goodIdx),y(goodIdx),'k.','MarkerSize',10)
%     xlabel('Ldir')
%     xlim([0 1])
%     ylabel('rPref')
% 
%     if svePlt == 1
%         figFileDir = fullfile(figDir,animal,exptName);
%         if ~isfolder(figFileDir)
%             mkdir(figFileDir)
%         end
%         figFileName = fullfile(figFileDir,[exptName '_p' num2str(probe) '_SummaryPlot']);
%         saveas(gcf,figFileName)
%     end

    
end

%% Save



end
