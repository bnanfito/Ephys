% visTest

function [pVis] = visTest(sumStats)

    for u = 1:height(sumStats)

        cndInclude = limitConds(sumStats(u,:));

        FR = sumStats(u,:).fr;
        R = sumStats.response{u};
        if sum(strcmp(sumStats.paramKey{u},'contrast')) == 1
            R = R(:,:,sumStats.oriPref(u));
        end
        Rblank = sumStats.rBlank{u};
        
        baseFR = FR.base(:,cndInclude);    baseFR = baseFR(:); baseFR = baseFR(~isnan(baseFR));
        stimFR = FR.stim(:,cndInclude);    stimFR = stimFR(:); stimFR = stimFR(~isnan(stimFR));
        bcFR = FR.bc(:,cndInclude);        bcFR = bcFR(:); bcFR = bcFR(~isnan(bcFR));
     
        pVis(u).signrank = signrank(bcFR);
        pVis(u).ranksum = ranksum(baseFR,stimFR);
        if ismember('rBlank',sumStats.Properties.VariableNames)
            pVis(u).anova = anova1([R,Rblank],[],'off');
        end

    end
    pVis = pVis';

end