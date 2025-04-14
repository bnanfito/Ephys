% visTest

function [pVis] = visTest(sumStats)

    for u = 1:height(sumStats)
    
        R = sumStats.response{u};

        % if stimulus size is a variable being changes, only consider
        % conditions that were used to calculate sumStats (ie the ones that
        % were used to make the response matrix, R)
        sizeId = contains(sumStats.paramKey{u},'size');
        cndSize = sumStats.cndKey{u}(:,sizeId);
        sizes = unique(cndSize);
        if sum(sizeId)>0 && length(sizes)>1
            ffR = sumStats(u,:).fr.bc(:,cndSize>100); %
            cmpre = R==ffR;
            cmpre = cmpre(~isnan(R) & ~isnan(ffR));
            isFF = length(unique(cmpre)) == 1 & max(unique(cmpre)) == 1;
            if isFF
                cndInclude = cndSize==max(sizes);
            else
                cndInclude = cndSize==min(sizes);
            end
        else
            cndInclude = ones(size(sumStats.cndKey{u},1),1)==1;
        end
        
        baseFR = sumStats(u,:).fr.base(:,cndInclude);    baseFR = baseFR(:); baseFR = baseFR(~isnan(baseFR));
        stimFR = sumStats(u,:).fr.stim(:,cndInclude);    stimFR = stimFR(:); stimFR = stimFR(~isnan(stimFR));
        bcFR = sumStats(u,:).fr.bc(:,cndInclude);        bcFR = bcFR(:); bcFR = bcFR(~isnan(bcFR));
     
        pVis(u).signrank = signrank(bcFR);
        pVis(u).ranksum = ranksum(baseFR,stimFR);
        pVis(u).anova = anova1([sumStats.response{u},sumStats.rBlank{u}],[],'off');

    end
    pVis = pVis';

end