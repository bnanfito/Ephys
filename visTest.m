% visTest

function [pVis] = visTest(sumStats)

    for u = 1:height(sumStats)

        cndInclude = limitConds(sumStats(u,:));
        
        baseFR = sumStats(u,:).fr.base(:,cndInclude);    baseFR = baseFR(:); baseFR = baseFR(~isnan(baseFR));
        stimFR = sumStats(u,:).fr.stim(:,cndInclude);    stimFR = stimFR(:); stimFR = stimFR(~isnan(stimFR));
        bcFR = sumStats(u,:).fr.bc(:,cndInclude);        bcFR = bcFR(:); bcFR = bcFR(~isnan(bcFR));
     
        pVis(u).signrank = signrank(bcFR);
        pVis(u).ranksum = ranksum(baseFR,stimFR);
        if ismember('rBlank',sumStats.Properties.VariableNames)
            pVis(u).anova = anova1([sumStats.response{u},sumStats.rBlank{u}],[],'off');
        end

    end
    pVis = pVis';

end