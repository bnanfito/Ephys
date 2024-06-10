% anaDotCoh

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

animal = 'febj5';
unit = '001';
expt = '006';
exptName = [animal '_u' unit '_' expt];
probe = 'PSS';

plt = 1;
plr = 1;
anaMode = 'SU';
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

cohInd = strcmp(trialInfo.dom,'dotCoherence');
cohVals = unique(trialInfo.domval(:,cohInd));

oriInd = strcmp(trialInfo.dom,'ori');
oriVals = unique(trialInfo.domval(:,oriInd));

xIdx = [fliplr(find(trialInfo.domval(:,oriInd)==oriVals(1))') find(trialInfo.domval(:,oriInd)==oriVals(2))'];

%% Compute 

[spks,trialExclude] = orgSpks(animal,unit,expt,probe,anaMode,dataFold);
[goodUnits,pVis] = screenUnits(spks,anaMode,blank,visTest,alpha);
nU = length(spks);

prefA = [];
prefB = [];
noPref = [];
for u = 1:nU

    exptID{u} = exptName;
    uID(u) = spks(u).unitId;
    probeID(u) = probe;
    areaID{u} = area;

    R{u} = spks(u).fr.bc(:,xIdx);
    C{u} = trialInfo.domval(xIdx,:)';
    Rb{u} = spks(u).fr.bc(:,blank);

    if ismember(u,find(goodUnits))
        
        if ttest(R{u}(:,1),R{u}(:,end))
            dMean = mean(R{u}(:,1),'omitnan') - mean(R{u}(:,end),'omitnan');
            if dMean>0
                prefA = [prefA;u];
            elseif dMean<0
                prefB = [prefB;u];
            end
        else
            noPref = [noPref;u];
        end

    end

end

for i = 1:length(xIdx)
    lbls{i} = [num2str(C{1}(2,i)) '%/' num2str(C{1}(1,i))];
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

        subplot(1,2,2);hold on;
        meanR = mean(R{u},'omitnan');
        sem = std(R{u},'omitnan')/sqrt(size(R{u},1));
        meanRb = mean(Rb{u},'omitnan');

        patch([1 length(xIdx) length(xIdx) 1],[min(Rb{u}) min(Rb{u}) max(Rb{u}) max(Rb{u})],'r','EdgeColor','none','FaceAlpha',0.2)
        yline(meanRb,'r--','LineWidth',1)
        plot(R{u}','.','Color',[0.8 0.8 0.8],'LineWidth',2)
        plot(repmat(1:length(xIdx),2,1),meanR+([-1;1]*sem),'k','LineWidth',2)
        plot(meanR,'k-o','LineWidth',2)
        xticks(1:length(xIdx))
        xticklabels(lbls)
        yline(0,'k')

        ttl = [exptName ' p' num2str(probe) ' (' area ') ' anaMode '#' num2str(uID(u))];
        if ~ismember(u,find(goodUnits))
            ttl = [ttl '(BAD UNIT)'];
        end
        sgtitle(ttl)
    end


    figure;hold on
    y = mean(mean(cat(3,R{prefA}),1,'omitnan'),3);
    plot(y,'r')
    y = mean(mean(cat(3,R{prefB}),1,'omitnan'),3);
    plot(y,'b')
    y = mean(mean(cat(3,R{noPref}),1,'omitnan'),3);
    plot(y,'k')
    
end

