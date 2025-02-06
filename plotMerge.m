
clear all
close all

%Settings
animalId = 'febj7';
mergeId = '000005000011000015';
probeId = 1;
mergeName = [animalId '_uMMM_' mergeId];
dataFold = 'F:\Brandon\data';
anaMode = 'SU';

plr = 0;
alignBit = 0;
visTest = 'ranksum';
alpha = 0.01;

%Load merge info
load(fullfile(dataFold,'Ephys',animalId,mergeName,[mergeName '_mergeInfo.mat']))
nFiles = length(mergeInfo.files);
for f = 1:nFiles
    exptName{f,1} = [animalId '_' mergeInfo.files{f}];

    %Load expt info
    load(fullfile(dataFold,'Ephys',animalId,exptName{f},[exptName{f} '_trialInfo.mat']))
    load(fullfile(dataFold,'Ephys',animalId,exptName{f},[exptName{f} '_id.mat']))
    load(fullfile(dataFold,'Ephys',animalId,exptName{f},[exptName{f} '.analyzer']),'-mat')

    varInfo = whos('probeId');
    if strcmp(varInfo.class,'char')
        area = probeId;
        probeId = find(strcmp({id.probes(:).area}',area));
    end
    clear varInfo
    area = id.probes(probeId).area;
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
        blankTrial = [];
    else
        blank = (1:nConds)==trialInfo.blankId;
        blankTrial = find(trialInfo.triallist==find(blank));
    end
    [sortTrialInd(:,1),sortTrialInd(:,2)] = sort(trialInfo.triallist);
    nT = ((1:(trialL*sf))-(predelay*sf))/sf;
    kernel = normpdf(-3:6/2000:3);
    colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
    
    oriInd = strcmp(trialInfo.dom,'ori');
    sizeInd = contains(trialInfo.dom,'size');
    posInd = strcmp(trialInfo.dom,'x_pos');
    if nDom==1 && sum(oriInd)==1 && sum(sizeInd)==0
        cndInclude = ~blank;
    elseif nDom == 2 && sum(sizeInd)==1 
        sizes = unique(trialInfo.domval(:,sizeInd));
        ffIdx = sizes>100;
        hemiIdx =  sizes<100;
        cndInclude = trialInfo.domval(:,sizeInd) == sizes(hemiIdx);
    elseif nDom == 2 && sum(posInd)==1 
        pos = unique(trialInfo.domval(:,posInd));
        ffIdx = pos==min(pos);
        hemiIdx = pos==max(pos);
        cndInclude = trialInfo.domval(:,sizeInd) == pos(hemiIdx);
    end
    trialInclude = ismember(trialInfo.triallist,find(cndInclude));

    %Compute spk info
    disp(['orgSpks: ' exptName{f}])
    [spks{f},~] = orgSpks(animalId,exptName{f}(8:10),exptName{f}(12:14),[num2str(probeId) '_' num2str(f)],'SU',dataFold);
    nU = length(spks{f});
    for u = 1:nU % u indexes a unit (column) in structure spks

        exptID{u} = exptName{f};
        uID(u) = spks{f}(u).unitId;
        probeID(u) = probeId;
        areaID{u} = area;
        paramKey{u} = trialInfo.dom;
        cndKey{u} = trialInfo.domval;
    
        if ~isempty(trialInfo.blankId)
            Rblank{u} = spks{f}(u).fr.stim(:,blank);
        end
    
        R{u}(:,:) = spks{f}(u).fr.bc(:,cndInclude);
        C{u}(:,:) = trialInfo.domval(cndInclude,:)';
    
        dir = C{u}(oriInd,:);
    
        ori = mod(dir,180);
        oris = unique(ori);
        for o = 1:length(oris)
            idx = ori == oris(o);
            rMean_ori(o) = mean(R{u}(:,idx),'all','omitnan');
        end
        rPref_ori = max(rMean_ori);
        if sum(rMean_ori==rPref_ori)>1
            rIn = rMean_ori;
            pks = rIn == rPref_ori;
            rConv = rIn([end 1:end 1]);
            rConv = conv(rConv,ones(1,3)*(1/3),'same');
            rConv = rConv(1+1:end-1);
            rConv(~pks) = 0;
            cPref_ori = oris(find(rConv==max(rConv),1,'first'));
            clear rIn pks rConv 
        else
            cPref_ori = oris(rMean_ori==rPref_ori);
        end
        cNull_ori = mod(cPref_ori+90,180);
        rNull_ori = rMean_ori(oris==cNull_ori);
        osi(u,1) = abs(rPref_ori-rNull_ori)/rPref_ori;
        clear rMean_ori rPref_ori rNull_ori
    
        rMean = mean(R{u}(:,:),'omitnan');
        rPref(u) = max(rMean);
        if sum(rMean==rPref(u))>1
            rIn = rMean;
            pks = rIn == rPref(u);
            rConv = rIn([end 1:end 1]);
            rConv = conv(rConv,ones(1,3)*(1/3),'same');
            rConv = rConv(1+1:end-1);
            rConv(~pks) = 0;
            cPref(u) = dir(find(rConv==max(rConv),1,'first'));
            clear rIn pks rConv 
        else
            cPref(u) = dir(rMean==rPref(u));
        end
        cNull(u) = mod(cPref(u)+180,360);
        rNull(u) = rMean(dir==cNull(u));
    
        dsi(u,1) = abs(rPref(u)-rNull(u))/rPref(u);
        mv{u} = meanvec(dir,rMean);
        ldr(u,1) = mv{u}.ldr;
        lor(u,1) = mv{u}.lor;
    
        %double gaussian fit
        G{u} = dirGauss(rMean,dir,0);
    
    end
    
    varNames = {'exptName','probe','area','uInfo','uID','spkTimes','latency','fr','paramKey','cndKey','response','condition','gaussFit','meanVec','rPref','oriPref','rNull','oriNull','rBlank','dsi','osi','ldr','lor'};
    sumStats{f} = table(exptID',probeID',areaID',{spks{f}.info}',uID',{spks{f}.stimCent}',vertcat(spks{f}.late),vertcat(spks{f}.fr),paramKey',cndKey',R',C',G',mv',rPref',cPref',rNull',cNull',Rblank',dsi,osi,ldr,lor,'VariableNames',varNames);
    
    [goodUnit,pVis] = screenUnits(sumStats{f},anaMode,blank,visTest,alpha);
    sumStats{f}.goodUnit = goodUnit';
    sumStats{f}.pVis = pVis;

    clear id trialInfo analyzer trialInclude
end
