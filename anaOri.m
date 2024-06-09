% anaOri

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
unit = '000';
expt = '023';
exptName = [animal '_u' unit '_' expt];
probe = 1;

plt = 1;
plr = 0;
anaMode = 'MU';
stimPres = 'ff';
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
else
    blank = (1:nConds)==trialInfo.blankId;
end
[~,sortTrialInd] = sort(trialInfo.triallist);
nT = ((1:(trialL*sf))-(predelay*sf))/sf;
kernel = normpdf(-3:6/2000:3);
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};


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

    if nDom==1 && strcmp(trialInfo.dom{1},'ori')

        cndInclude = ~blank;

    elseif nDom==2 && sum(strcmp(trialInfo.dom,'y_size'))==1

        if strcmp(stimPres,'ff')
            cndInclude = trialInfo.domval( : , strcmp(trialInfo.dom,'y_size') ) >=150;
        elseif strcmp(stimPres,'hemi')
            cndInclude = trialInfo.domval( : , strcmp(trialInfo.dom,'y_size') ) <=75;
        end

    end

    if ~isempty(trialInfo.blankId)
        Rblank{u} = spks(u).fr.bc(:,blank);
    end
    R{u} = spks(u).fr.bc(:,cndInclude);
    C{u} = trialInfo.domval(cndInclude, strcmp(trialInfo.dom,'ori') )';
    rMean = mean(R{u},'omitnan');
    rPref = max(rMean);
    if sum(rMean==rPref)>1
        rIn = rMean;
        pks = rIn == rPref;
        rConv = rIn([end 1:end 1]);
        rConv = conv(rConv,ones(1,3)*(1/3),'same');
        rConv = rConv(1+1:end-1);
        rConv(~pks) = 0;
        cPref = C{u}(find(rConv==max(rConv),1,'first'));
        clear rIn pks rConv 
    else
        cPref = C{u}(rMean==rPref);
    end
    cNull = mod(cPref+180,360);
    rNull = rMean(C{u}==cNull);

    dsi(u) = abs(rPref-rNull)/rPref;
    mv = meanvec(C{u},rMean);
    ldr(u) = mv.l;

    Rpref(u) = rPref;
    Cpref(u) = cPref;
    Rnull(u) = rNull;
    Cnull(u) = cNull;

end

varNames = {'exptName','probe','area','uID','goodUnit','pVis','fr','response','condition','rPref','cPref','rNull','cNull','rBlank','dsi','ldr'};
sumStats = table(exptID',probeID',areaID',uID',goodUnits',pVis,vertcat(spks.fr),R',C',Rpref',Cpref',Rnull',Cnull',Rblank',dsi',ldr','VariableNames',varNames);

if strcmp(anaMode,'MU')
    sumStats.xPos = vertcat(spks.xPos);
    sumStats.zPos = vertcat(spks.zPos);
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

        if plr == 1
            subplot(1,2,2,polaraxes);hold on
            x = deg2rad(C{u}); x = [x(1:end) x(1)];
            y = mean(R{u},'omitnan'); y = [y(1:end) y(1)];
            polarplot(x,y,'k-o')            
        else
            subplot(1,2,2);hold on
            plot(C{u},mean(R{u},'omitnan'),'k-o')
        end

        ttl = [exptName ' p' num2str(probe) ' (' area ') ' anaMode '#' num2str(u)];
        if ~ismember(u,find(goodUnits))
            ttl = [ttl '(BAD UNIT)'];
        end
        sgtitle(ttl)
    end

    
end

%% Save


