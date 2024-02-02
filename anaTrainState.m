%anaStateSpace

clear all
close all

if ispc
%     dataFold = 'D:\data'; 
    dataFold = 'F:\Brandon\data';
elseif ismac
    dataFold = '/Volumes/Lab drive/Brandon/data';
%     dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end
physDir = fullfile(dataFold,'Ephys');

area = 'PSS';

% % V1 COOLED
% animals{1} = 'febh2';
% animals{2} = 'febh3';
% animals{3} = 'febj3';

% CONTROL (AUGUSTO)
animals{1} = 'FEAO4';
animals{2} = 'FEAQ5';
animals{3} = 'FEAS6';
animals{4} = 'FEAT1';
animals{5} = 'FEAN6';

nAnimals = length(animals);
B = [];
B_condID = [];
A = [];
A_condID = [];
for a = 1:nAnimals

    if strcmp(animals{a},'febh2')
        before = {'_u000_001'};
        after = {'_u000_040','_u000_041','_u002_000'};
        trainAx = [90 270];
    elseif strcmp(animals{a},'febh3')
        before = {'_u000_002'};
        after = {'_u000_043','_u001_000','_u002_002','_u003_000'};
        trainAx = [0 180];
    elseif strcmp(animals{a},'febj3')
        before = {'_u000_001'};
        after = {'_u000_038'};
        trainAx = [0 180];
    elseif strcmp(animals{a},'febj2')
        before = {'_u000_003'};
        after = {'_u000_023'};
        trainAx = [90 270];
    elseif strcmp(animals{a},'febj4')
        before = {'_u000_006'};
        after = {'_u000_027'};
        trainAx = [90 270];
    elseif strcmp(animals{a},'FEAO4')
        before = {'_u000_001'};
        after = {'_u002_000'};
        trainAx = [0 180];
    elseif strcmp(animals{a},'FEAQ5')
        before = {'_u000_000'};
        after = {'_u002_000'};
        trainAx = [90 270];
    elseif strcmp(animals{a},'FEAS6')
        before = {'_u000_003'};
        after = {'_u002_001'};
        trainAx = [60 240];
    elseif strcmp(animals{a},'FEAT1')
        before = {'_u000_000'};
        after = {'_u002_000'};
        trainAx = [90 270];
    elseif strcmp(animals{a},'FEAN6')
        before = {'_u000_000'};
        after = {'_u002_002'};
        trainAx = [150 330];
    elseif strcmp(animals{a},'FEAQ2')
        before = {'_u000_001'};
        after = {'_u002_001'};
        trainAx = [];
    elseif strcmp(animals{a},'FEAQ3')
        before = {'_u000_001'};
        after = {'_u001_001'};
        trainAx = [];
    elseif strcmp(animals{a},'FEAQ4')
        before = {'_u000_001'};
        after = {'_u001_001'};
        trainAx = [];
    elseif strcmp(animals{a},'FEAQ7')
        before = {'_u000_000'};
        after = {'_u002_001'};
        trainAx = [];
    end
    
    exptNameB = [animals{a} before{1}];unit = exptNameB(8:10); expt = exptNameB(12:14);
    load(fullfile(physDir,animals{a},exptNameB,[exptNameB '_id.mat']))
    for p = 1:length(id.probes)
        if strcmp(id.probes(p).area,area)
            [b,b_condID] = stateSpace(animals{a},unit,expt,p,0);
        end
    end
    B = horzcat(B,b);
    B_condID = horzcat(B_condID,b_condID);

    exptNameA = [animals{a} after{1}];unit = exptNameA(8:10); expt = exptNameA(12:14);
    load(fullfile(physDir,animals{a},exptNameA,[exptNameA '_id.mat']))
    for p = 1:length(id.probes)
        if strcmp(id.probes(p).area,area)
            [a,a_condID] = stateSpace(animals{a},unit,expt,p,0);
        end
    end
    A = horzcat(A,a);
    A_condID = horzcat(A_condID,a_condID);

end

figure;
for tr = 1:2

    if tr==1
        X = B;
        X_condID = mean(B_condID,2);
    elseif tr==2
        X = A;
        X_condID = mean(A_condID,2);
    end
    
    [COEFF, SCORE] = pca(X);
    
    subplot(1,2,tr); hold on
    oris = [0 90 180 270];
    % colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840],[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840],[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
    colors = {'r','g','b',[0.8500 0.3250 0.0980]};
    countC = 0;
%     if isempty(trialInfo.blankId)
%         conds = trialInfo.domval;
%     else
        conds = [unique(X_condID);-1];
%     end
    for c = 1:length(conds)
        cond = conds(c);
        if ~ismember(cond,oris)
            continue
        end
        countC = countC+1;
        plot3(SCORE(X_condID==cond,1),SCORE(X_condID==cond,2),SCORE(X_condID==cond,3),'-o','Color',colors{countC})
    
    end

end