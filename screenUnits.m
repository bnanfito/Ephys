function [out,pVis] = screenUnits(sumStats,mode,blankBit,visTest,alpha)

% INPUTS
% 1. in = structure containing spike data organized by unit.
% 2. mode = 'SU' for single-unit analysis or 'MU' for multi-unit analysis
% 3. blankBit = binary vector (0 or 1) indicating which condition in the nRep x nCondition firing rate matrix (field in 'in') is a blank (=1)

% OUTPUTS
% out = nUnits long binary vector (0 or 1) indicating which pass inclusion criteria (=1) and which do not (=0)
% pVis = vector of p values for each unit, calculated from hypothesis test (visTest)
 
% alpha = 0.01;
% visTest = 'ranksum';

for u = 1:height(sumStats)

    if length(sumStats.paramKey{u})>1
        sizeInd = contains(sumStats.paramKey{u},'size');
        sizeVec = sumStats.cndKey{u}(:,sizeInd);
        cndInclude = find(sizeVec==min(sizeVec));
    else
        cndInclude = ~blankBit;
    end

    baseFR = sumStats(u,:).fr.base(:,cndInclude);    baseFR = baseFR(:); 
    stimFR = sumStats(u,:).fr.stim(:,cndInclude);    stimFR = stimFR(:); 
    bcFR = sumStats(u,:).fr.bc(:,cndInclude);        bcFR = bcFR(:); bcFR = bcFR(~isnan(bcFR));
 
    if strcmp(visTest,'signrank')
        pVis(u,1) = signrank(bcFR);
    elseif strcmp(visTest,'ranksum')
        pVis(u,1) = ranksum(baseFR,stimFR);
    elseif strcmp(visTest,'anova')
        pVis(u,1) = anova1(sumStats(u,:).fr.bc,[],'off');
    end
    isVis(u) = pVis(u,1)<alpha;

    if strcmp(mode,'SU')
        isAct(u) = max(mean(sumStats(u,:).fr.bc(:,cndInclude),'omitnan'))>=2;
        isSU(u) = strcmp(sumStats(u,:).uInfo,'SU');
        maybeSU(u) = strcmp(sumStats(u,:).uInfo,'SU?');
    elseif strcmp(mode,'MU')
        isAct(u) = max(mean(sumStats(u,:).fr.bc(:,cndInclude),'omitnan'))>=2;
    end

end

if strcmp(mode,'SU')
    out = isVis & isAct & isSU;
    %     out = isVis & isAct & (isSU | maybeSU);
elseif strcmp(mode,'MU')
    out = isVis & isAct;
end



end