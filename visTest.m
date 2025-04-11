% visTest

function [pVis] = visTest(sumStats)

    for u = 1:height(sumStats)
    
        if sum(contains(sumStats.paramKey{u},'size'))>0
            cndSize = sumStats.cndKey{u}(:,contains(sumStats.paramKey{u},'size'));
            sizes = unique(cndSize);
            if length(sizes)>1
                a = sumStats.response{u};
                b = sumStats(u,:).fr.bc(:,cndSize>100);
                cmpre = a==b;
                cmpre = cmpre(~isnan(a) & ~isnan(b));
                isFF = length(unique(cmpre)) == 1 & max(unique(cmpre)) == 1;
                if isFF
                    cndInclude = cndSize==max(sizes);
                else
                    cndInclude = cndSize==min(sizes);
                end
            else
                cndInclude = ~blankBit;
            end
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
            pVis(u,1) = anova1([sumStats.response{u},sumStats.rBlank{u}],[],'off');
        end
        isVis(u) = pVis(u,1)<alpha;
    
        if strcmp(mode,'SU')
            isAct(u) = sumStats.rPref(u)>=2;
            isSU(u) = strcmp(sumStats(u,:).uInfo,'SU');
            maybeSU(u) = strcmp(sumStats(u,:).uInfo,'SU?');
        elseif strcmp(mode,'MU')
            isAct(u) = sumStats.rPref(u)>=2;
        end
    
    end

end