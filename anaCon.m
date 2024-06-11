% anaContrast

%% Initialize

clear
close all

if ispc
    dataFold = 'D:\data'; 
%     dataFold = 'C:\Users\brand\Documents\data';
%     dataFold = 'F:\Brandon\data';
elseif ismac
    dataFold = '/Volumes/Lab drive/Brandon/data';
%     dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end
physDir = fullfile(dataFold,'Ephys');
figDir = fullfile(dataFold,'Figures');


%% Settings

animal = 'febl0';
unit = '000';
expt = '003';
exptName = [animal '_u' unit '_' expt];
probe = 'V1';

plt = 1;
plr = 1;
anaMode = 'MU';
stimPres = 'ff';
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
else
    blank = (1:nConds)==trialInfo.blankId;
end
[~,sortTrialInd] = sort(trialInfo.triallist);
nT = ((1:(trialL*sf))-(predelay*sf))/sf;
kernel = normpdf(-3:6/2000:3);
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};

conInd = strcmp(trialInfo.dom,'contrast');
conVals = unique(trialInfo.domval(:,conInd));

oriInd = strcmp(trialInfo.dom,'ori');
oriVals = unique(trialInfo.domval(:,oriInd));


%% Compute 

[spks,trialExclude] = orgSpks(animal,unit,expt,probe,anaMode,dataFold);
[goodUnits,pVis] = screenUnits(spks,anaMode,blank,visTest,alpha);
nU = length(spks);

for u = 1:nU

    exptID{u} = exptName;
    uID(u) = spks(u).unitId;
    probeID(u) = probe;
    areaID{u} = area;

    if nDom==2 && sum(strcmp(trialInfo.dom,'ori'))==1 && sum(strcmp(trialInfo.dom,'contrast'))==1

        cndInclude = ~blank;

    elseif nDom==3 && sum(strcmp(trialInfo.dom,'y_size '))==1

        if strcmp(stimPres,'ff')
            cndInclude = trialInfo.domval( : , strcmp(trialInfo.dom,'y_size ') ) >=150;
        elseif strcmp(stimPres,'hemi')
            cndInclude = trialInfo.domval( : , strcmp(trialInfo.dom,'y_size ') ) <=75;
        end

    end

    if ~isempty(trialInfo.blankId)
        Rb{u} = spks(u).fr.bc(:,blank);
    end
    for o = 1:length(oriVals)

        curOri = oriVals(o);
        curOriIdx = (trialInfo.domval(:,oriInd) == curOri)';
        if ~isempty(trialInfo.blankId)
            curOriIdx = [curOriIdx false];
        end
        curCondInclude = cndInclude & curOriIdx;

        R{u,o} = spks(u).fr.bc(:,curCondInclude);
        C{u,o} = trialInfo.domval(curCondInclude, conInd )';

    end
    
end

for i = 1:length(C{1,1})
    lbls{i} = num2str(C{1,1}(i));
end



%% Plot

if plt == 1

    for u = 1:nU
        figure;

        subplot(2,2,1);hold on
        bins = -1:0.1:2;
        h = histogram(spks(u).stimCent(1,:),'BinEdges',bins);
        h.FaceColor = 'k';
        h.EdgeColor = 'none';
        xlim([bins(1) bins(end)])

        subplot(2,2,3);hold on
        x = spks(u).stimCent(1,:);
        y = spks(u).stimCent(2,:);
        patch([0 1 1 0],[0 0 max(y)+1 max(y)+1],'k','EdgeColor','none','FaceAlpha',0.2)
        plot(x,y,'k.')
        xlim([-1 2])
        ylim([0 max(y)+1])

        subplot(1,2,2);hold on
        patch([1 length(C{1,1}) length(C{1,1}) 1],[min(Rb{u}) min(Rb{u}) max(Rb{u}) max(Rb{u})],'r','EdgeColor','none','FaceAlpha',0.2)
        meanRb = mean(Rb{u},'omitnan');
        yline(meanRb,'r--','LineWidth',1)
        for o = 1:length(oriVals)
            meanR = mean(R{u,o},'omitnan');
            sem = std(R{u,o},'omitnan')/sqrt(size(R{u},1));

            plot(R{u,o}','.','Color',[0.8 0.8 0.8],'LineWidth',2)
            plot(repmat(1:length(C{1,1}),2,1),meanR+([-1;1]*sem),'k','LineWidth',2)
            plot(meanR,'k-o','LineWidth',2)
        end
        xticks(1:length(C{1,1}))
        xticklabels(lbls)
        yline(0,'k')


        ttl = [exptName ' p' num2str(probe) ' (' area ') ' anaMode '#' num2str(uID(u))];
        if ~ismember(u,find(goodUnits))
            ttl = [ttl '(BAD UNIT)'];
        end
        sgtitle(ttl)
    end


    
end


