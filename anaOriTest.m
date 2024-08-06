
%test - anaOri

function [sumStats] = anaOriTest(animal,unit,expt,probe,anaMode,dataFold,plt,svePlt)

    %% settings
    plr = 0;
    alignBit = 1;
    stimMode = 'mono c hemi';
    visTest = 'ranksum';
    alpha = 0.01;
    
    %% load 
    exptName = [animal '_u' unit '_' expt];
    physDir = fullfile(dataFold,'Ephys');
    figDir = fullfile(dataFold,'Figures');
    exptDir = fullfile(physDir,animal,exptName);
    load(fullfile(exptDir, [exptName '_id.mat'] ))
    load(fullfile(exptDir, [exptName '_trialInfo.mat'] ))
    load(fullfile(exptDir, [exptName '.analyzer']),'-mat')

    %% initialize
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
    if stimTime < 1
        st = stimTime;
    else
        st = 1;
    end
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
    
    [sortTrialInd(:,1),sortTrialInd(:,2)] = sort(trialInfo.triallist);
    nT = ((1:(trialL*sf))-(predelay*sf))/sf;
    kernel = normpdf(-3:6/2000:3);
    colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]...
             ,[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330]...
             ,[0.6350 0.0780 0.1840]};

    if strcmp(Analyzer.modID,'PG')
        stimStartID = 3;
    elseif strcmp(Analyzer.modID,'RD')
        stimStartID = 7;
    end

    if ismember('contrast',trialInfo.dom)
        xLab = 'contrast';
    elseif ismember('s_freq',trialInfo.dom)
        xLab = 'spatial freq.';
    else
        xLab = 'orientation';
    end

    if ismember('y_size',trialInfo.dom) || ismember('x_size',trialInfo.dom) || ismember('y_size ',trialInfo.dom)
        if strcmp(stimMode,'mono c hemi')
            cndInclude = trialInfo.dom
        elseif strcmp(stimMode,'bi ff')

        end
    else
    
    end

    
    %% compute
    [spks] = orgSpks(animal,unit,expt,probe,anaMode,dataFold);
    [goodUnits,pVis] = screenUnits(spks,anaMode,blank,visTest,alpha);
    
    
    
    
    
    
    %% plotting
    if plt == 1
        figure; hold on
    end


end
