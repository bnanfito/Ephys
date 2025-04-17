function [out] = screenUnits(sumStats,anaMode,varargin)
 
alpha = 0.01;
% if a statistical test for visual responsiveness is not specified in the
% 3rd argin, ranksum will be used by default
if isempty(varargin)
    visTest = 'ranksum';
else
    visTest = varargin{1};
end

for u = 1:height(sumStats)

    switch visTest
        case 'anova'
            isVis = sumStats.pVis(u).anova<alpha;
        case 'ranksum'
            isVis = sumStats.pVis(u).ranksum<alpha;
        case 'signrank'
            isVis = sumStats.pVis(u).signrank<alpha;
    end
    
    if strcmp(anaMode,'SU')
        isAct = sumStats.rPref(u)>=2;
        isSU = strcmp(sumStats.uInfo{u},'SU');
        maybeSU = strcmp(sumStats.uInfo{u},'SU?');
        out(u,1) = isVis & isAct & isSU;
%         out(u,1) = isVis & isAct & (isSU | maybeSU);
    elseif strcmp(anaMode,'MU')
        isAct = sumStats.rPref(u)>=2;
        out(u,1) = isVis & isAct;
    end

    clear isVis isAct isSU maybeSU
end



end