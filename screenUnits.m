function [out] = screenUnits(sumStats,anaMode)
 
alpha = 0.01;
for u = 1:height(sumStats)

    isVis = sumStats.pVis(u).anova<alpha;
    if strcmp(anaMode,'SU')
        isAct = sumStats.rPref(u)>=2;
        isSU = strcmp(sumStats{u}.uInfo,'SU');
        maybeSU = strcmp(sumStats{u}.uInfo,'SU?');
        out(u) = isVis & isAct & isSU;
%         out(u) = isVis & isAct & (isSU | maybeSU);
    elseif strcmp(anaMode,'MU')
        isAct = sumStats.rPref(u)>=2;
        out(u) = isVis & isAct;
    end

    clear isVis isAct isSU maybeSU
end



end