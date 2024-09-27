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
visTest = 'anova';
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
posInd = strcmp(trialInfo.dom,'x_pos');
if nDom==1 && strcmp(trialInfo.dom{1},'ori') && sum(sizeInd)==0
    nStim = 1;
    stTrialList{1} = find(trialInfo.triallist~=find(blank));
elseif nDom == 2 && sum(sizeInd)>0 
    sizes = unique(trialInfo.domval(:,sizeInd));
    nStim = length(sizes);
    for st = 1:nStim
        stTrialList{st} = find(ismember(trialInfo.triallist,find(trialInfo.domval(:,sizeInd)==sizes(st))));
    end
elseif nDom == 2 && sum(posInd)>0 
    pos = unique(trialInfo.domval(:,posInd));
    nStim = length(pos);
    for st = 1:nStim
        stTrialList{st} = find(ismember(trialInfo.triallist,find(trialInfo.domval(:,posInd)==pos(st))));
    end
end

%% Compute 

[spks,trialExclude] = orgSpks(animal,unit,expt,probe,anaMode,dataFold);
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

    for st = 1:nStim

        if nDom==1 && strcmp(trialInfo.dom{1},'ori') && sum(sizeInd)==0
    
            cndInclude = ~blank;
    
        elseif nDom == 2 && sum(sizeInd)>0
            cndInclude = trialInfo.domval( : , sizeInd ) == sizes(st);
            if sizes(st) >= 150
                ffIdx = st;
                continue
            elseif sizes(st) <= 75
                hemiIdx = st;
            end
        elseif nDom == 2 && sum(posInd)>0
            cndInclude = trialInfo.domval( : , posInd ) == pos(st);
%             if pos(st) >= 
%                 ff = st;
%             elseif pos(st) <= 
%                 hemi = st;
%             end
        end
    
        if ~isempty(trialInfo.blankId)
            Rblank{u} = spks(u).fr.stim(:,blank);
        end
        R{u}(:,:,st) = spks(u).fr.bc(:,cndInclude);
        C{u}(:,:,st) = trialInfo.domval(cndInclude,:)';
        dir = C{u}(oriInd,:,st);
        paramKey{u} = trialInfo.dom;
        cndKey{u} = trialInfo.domval;
        
        rMean = mean(R{u}(:,:,st),'omitnan');
        rPref(u,1) = max(rMean);
        if sum(rMean==rPref(u,1))>1
            rIn = rMean;
            pks = rIn == rPref(u,1);
            rConv = rIn([end 1:end 1]);
            rConv = conv(rConv,ones(1,3)*(1/3),'same');
            rConv = rConv(1+1:end-1);
            rConv(~pks) = 0;
            cPref(u,1) = dir(find(rConv==max(rConv),1,'first'));
            clear rIn pks rConv 
        else
            cPref(u,1) = dir(rMean==rPref(u,1));
        end
        cNull(u,1) = mod(cPref(u,1)+180,360);
        rNull(u,1) = rMean(C{u}(oriInd,:,st)==cNull(u,1));
    
        dsi(u,1) = abs(rPref(u,1)-rNull(u,1))/rPref(u,1);
        mv{u} = meanvec(C{u}(oriInd,:,st),rMean);
        ldr(u,1) = mv{u}.ldr;

        %double gaussian fit
        G{u} = dirGauss(rMean,dir,0);


    end

end

varNames = {'exptName','probe','area','uInfo','uID','spkTimes','latency','fr','paramKey','cndKey','response','condition','gaussFit','meanVec','rPref','oriPref','rNull','oriNull','rBlank','dsi','ldr'};
sumStats = table(exptID',probeID',areaID',{spks.info}',uID',{spks.stimCent}',vertcat(spks.late),vertcat(spks.fr),paramKey',cndKey',R',C',G',mv',rPref,cPref,rNull,cNull,Rblank',dsi,ldr,'VariableNames',varNames);

[goodUnit,pVis] = screenUnits(sumStats,anaMode,blank,visTest,alpha);
sumStats.goodUnit = goodUnit';
sumStats.pVis = pVis;

if strcmp(anaMode,'MU')
    sumStats.xPos = vertcat(spks.xPos);
    sumStats.zPos = vertcat(spks.zPos);
end


%% Plot

if plt == 1

    for u = find(goodUnit)

        
        figure;
        legLbl{nStim+1} = 'blank';
        x = sumStats.spkTimes{1}(1,:);
        y = sumStats.spkTimes{1}(2,:);
        blankSpkIdx = ismember(y,blankTrialList);
        for st = 1:nStim

            clr = colors{st};
            if nDom==2 && sum(sizeInd)>0
                if st==ffIdx
                    legLbl{st} = ['full field'];
                elseif st==hemiIdx
                    legLbl{st} = ['hemi field'];
                end
            else
                legLbl{st} = ['data'];
            end

            szSpkIdx = ismember(y,stTrialList{st});

            subplot(2,2,1);hold on
            bins = -1:0.1:2;
            h = histogram(x(szSpkIdx),'BinEdges',bins);
            h.FaceColor = colors{st};
            h.EdgeColor = 'none';
            if ~isnan(sumStats.latency(u))
                xline(sumStats.latency(u),'--','LineWidth',2)
            end
            xlim([bins(1) bins(end)])
    
            subplot(2,2,3);hold on
            if isempty(x) || isempty(y)
                text(0,0,'no spikes')
                ttl = [exptName ' p' num2str(probe) ' (' area ') ' anaMode '#' num2str(uID(u))];
                if ~ismember(u,find(goodUnit))
                    ttl = [ttl '(BAD UNIT)'];
                end
                sgtitle(ttl)
                continue
            end
            if st == 1
                patch([0 1 1 0],[0 0 max(y)+1 max(y)+1],'k','EdgeColor','none','FaceAlpha',0.2)
            end
            for t = 1:nTrials
                if ismember(t,find(trialExclude))
                    patch([-predelay stimTime+postdelay stimTime+postdelay -predelay],[t-0.5 t-0.5 t+0.5 t+0.5],'r','EdgeColor','none','FaceAlpha',0.2)
                end
            end
            plot(x(szSpkIdx),y(szSpkIdx),'.','Color',colors{st})
            if ~isnan(sumStats.latency(u))
                xline(sumStats.latency(u),'--','LineWidth',2)
            end
            xlim([-1 2])
            ylim([0 max(y)+1])

            xP = sumStats.condition{u}(oriInd,:,st);
            yP = sumStats.response{u}(:,:,st);
            yMean = mean(yP,'omitnan');
            if alignBit == 1
                [xP,yMean,i] = alignDirTuning(xP,yMean);
                yP{u} = yP{u}(:,i);
            else
                xP = [xP xP(1)+360];
                yP = [yP yP(:,1)];
                yMean = [yMean yMean(1)];
            end
            sem = std(yP,'omitnan')/sqrt(size(yP,1));

            if plr == 1
                subplot(1,2,2,polaraxes);hold on
                xP = deg2rad(xP);
                polarplot(xP,yMean,'o','Color',clr);
                hold on;
%                 polarplot(xP,yP,'.','Color',clr)
                polarplot(repmat(xP,2,1),yMean+([1;-1]*sem),'Color',clr)
                
                polarplot(repmat(deg2rad(sumStats.meanVec{u}.angDir),2,1),...
                    [0 sumStats.meanVec{u}.magDir],'k','LineWidth',2)
                polarplot(repmat(deg2rad(sumStats.meanVec{u}.angDir),2,1),...
                    [0 sumStats.meanVec{u}.ldr*sumStats.meanVec{u}.magDir],'g','LineWidth',2)
                polarplot(deg2rad(sumStats.oriPref(u)),sumStats.rPref(u),'r*')
            else
                subplot(1,2,2);hold on

                plot(xP,yP','.','Color',clr)
                plot(repmat(xP,2,1),yMean+([1;-1]*sem),'Color',clr)
                plot(xP,yMean,'o','Color',clr);

                plot(repmat(sumStats.meanVec{u}.angDir,2,1),[0;sumStats.meanVec{u}.magDir],'k','LineWidth',2)
                plot(repmat(sumStats.meanVec{u}.angDir,2,1),[0;sumStats.meanVec{u}.ldr*sumStats.meanVec{u}.magDir],'g','LineWidth',2)
                plot(sumStats.oriPref(u),sumStats.rPref(u),'r*')
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

        ttl = [exptName ' p' num2str(probe) ' (' area ') ' anaMode '#' num2str(uID(u))];
        if ~ismember(u,find(goodUnit))
            ttl = [ttl '(BAD UNIT)'];
        end
        sgtitle(ttl)
        if svePlt == 1 && ismember(u,find(goodUnit))
            figFileDir = fullfile(figDir,animal,exptName);
            if ~isfolder(figFileDir)
                mkdir(figFileDir)
            end
            figFileName = fullfile(figFileDir,[exptName '_p' num2str(probe) '_' anaMode num2str(uID(u))]);
            saveas(gcf,figFileName)
        end
    end

    
end



end
