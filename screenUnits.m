function [out,pVis] = screenUnits(in,mode,blankBit,visTest,alpha)

% INPUTS
% 1. in = structure containing spike data organized by unit. Output of 'orgSpks' function
% 2. mode = 'SU' for single-unit analysis or 'MU' for multi-unit analysis
% 3. blankBit = binary vector (0 or 1) indicating which condition in the nRep x nCondition firing rate matrix (field in 'in') is a blank (=1)

% OUTPUTS
% out = nUnits long binary vector (0 or 1) indicating which pass inclusion criteria (=1) and which do not (=0)
% pVis = vector of p values for each unit, calculated from hypothesis test (visTest)
 
% alpha = 0.01;
% visTest = 'ranksum';

for u = 1:length(in)

    baseFR = in(u).fr.base(:,~blankBit);    baseFR = baseFR(:); 
    stimFR = in(u).fr.stim(:,~blankBit);    stimFR = stimFR(:); 
    bcFR = in(u).fr.bc(:,~blankBit);        bcFR = bcFR(:); bcFR = bcFR(~isnan(bcFR));
 
    if strcmp(visTest,'signrank')
        pVis(u,1) = signrank(bcFR);
    elseif strcmp(visTest,'ranksum')
        pVis(u,1) = ranksum(baseFR,stimFR);
    elseif strcmp(visTest,'anova')
        pVis(u,1) = anova1(in(u).fr.bc,[],'off');
    end
    isVis(u) = pVis(u,1)<alpha;

    if strcmp(mode,'SU')
        isAct(u) = max(mean(in(u).fr.bc(:,~blankBit),'omitnan'))>=2;
        isSU(u) = strcmp(in(u).info,'SU');
        maybeSU(u) = strcmp(in(u).info,'SU?');
    elseif strcmp(mode,'MU')
        isAct(u) = max(mean(in(u).fr.bc(:,~blankBit),'omitnan'))>=2;
    end

end

if strcmp(mode,'SU')
    out = isVis & isAct & isSU;
    %     out = isVis & isAct & (isSU | maybeSU);
elseif strcmp(mode,'MU')
    out = isVis & isAct;
end



end