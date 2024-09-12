clear all
close all

animal = 'febn1';
unit = {'001','001','001'};
expt = {'000','001','002'};
coolCnd = [1,2,3];
area = 'PSS';
anaMode = 'MU';
stimMode = 'hemi';
dataFold = '/Volumes/Lab drive/Brandon/data';

for e = 1:length(expt)
    [sumStats{e}] = anaOri(animal,unit{e},expt{e},area,anaMode,dataFold,0,0);
end

figure; hold on
for e = 1:length(sumStats)

    tbl = sumStats{e};

    sizeIdx = contains(tbl.paramKey{1},'size');
    if sum(sizeIdx) > 0
        sizes = unique(tbl.cndKey{1}(:,sizeIdx));
        if strcmp(stimMode,'hemi')
            stimIdx = sizes<=75;
        else strcmp(stimMode,'ff')
            stimIdx = sizes>=150;
        end
        cndIdx = tbl.cndKey{1}(:,sizeIdx) == sizes(stimIdx);
    else

    end


    if coolCnd(e) == 1
        clr = 'k';
    elseif coolCnd(e) == 2
        clr = 'c';
    elseif coolCnd(e) == 3
        clr = 'm';
    end
    
    isAct = tbl.rPref(:,stimIdx)>2;
    for u = 1:height(tbl)
%         pVis(u) = anova1([tbl.response{u} tbl.rBlank{u}],[],"off");
        baseFR = tbl.fr(u).base(:,cndIdx);
        stimFR = tbl.fr(u).stim(:,cndIdx);
        pVis(u) = ranksum(baseFR(:),stimFR(:));
        clear baseFR stimFR
    end
    isVis = pVis<0.01; clear pVis
    goodIdx = isAct' & isVis;

    cdf = cdfplot(tbl.latency(goodIdx));
    cdf.Color = clr;

end
