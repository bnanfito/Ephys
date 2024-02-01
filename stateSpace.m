%state space

function [X,X_condID] = stateSpace(animal,unit,expt,probe)

if ispc
%     dataFold = 'D:\data'; 
    dataFold = 'F:\Brandon\data';
elseif ismac
    dataFold = '/Volumes/Lab drive/Brandon/data';
%     dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end
physDir = fullfile(dataFold,'Ephys');
% animal = 'FEAN6';
% unit = '000';
% expt = '000';
% probe = 1;

exptName = [animal '_u' unit '_' expt];
spks = orgSpks(animal,unit,expt,probe,'MU',dataFold);
load(fullfile(physDir,animal,exptName,[exptName '_trialInfo.mat']))

for u = 1:length(spks)

    h = spks(u).psth_stim.values;
    H = [];
    condID = [];
    for c = 1:size(h,3)
        hC = sum(h(:,:,c));
        nReps = size(hC,1);
        nBins = size(hC,2);
        if c == trialInfo.blankId
            condID = vertcat(condID,repmat(-1,1,nBins));
        else
            condID = vertcat(condID,repmat(trialInfo.domval(c),1,nBins));
        end
        H = vertcat(H,hC);
    end

    popAct_condID(:,:,u) = condID;
    popAct(:,:,u) = H; % dim1 = nConds ; dim2  = nBins of psth
%     figure;imagesc(H);

    clear h H condID

end

obs = 0;
for t = 1:size(popAct,1)

    for b = 1:size(popAct,2)

        obs = obs+1;
        X_condID(obs,1) = unique(popAct_condID(t,b,:));
        X(obs,:) = popAct(t,b,:);

    end

end

% [COEFF, SCORE] = pca(X);
% 
% figure; hold on
% oris = [0 90 180 270];
% % colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840],[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840],[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
% colors = {'r','g','b','y'};
% countC = 0;
% if isempty(trialInfo.blankId)
%     conds = trialInfo.domval;
% else
%     conds = [trialInfo.domval;-1];
% end
% for c = 1:length(conds)
%     cond = conds(c);
%     if ~ismember(cond,oris)
%         continue
%     end
%     countC = countC+1;
%     plot3(SCORE(X_condID==c,1),SCORE(X_condID==c,2),SCORE(X_condID==c,3),'-o','Color',colors{countC})
% 
% end



end



