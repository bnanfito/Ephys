% plotSum

function [dist] = plotSum(sumStats,stimSize,plt)

    stimSizeInd = ismember(sumStats.paramKey{1,:},'x_size') | ...
                    ismember(sumStats.paramKey{1,:},'y_size') | ...
                    ismember(sumStats.paramKey{1,:},'y_size ');
    if sum(stimSizeInd)>0
        ffIdx = find(squeeze(sumStats.condition{1}(stimSizeInd,1,:))>100);
        hemiIdx = find(squeeze(sumStats.condition{1}(stimSizeInd,1,:))<100);
        if strcmp(stimSize,'hemi')
            setIdx = hemiIdx;
        elseif strcmp(stimSize,'ff')
            setIdx = ffIdx;
        end
    else
        setIdx = 1;
    end

    nU = height(sumStats);
    count = 0;
    for u = 1:nU

        R = sumStats.response{u}(:,:,setIdx);

        isAct = sumStats.rPref(u,setIdx)>2;
%         isVis = anova1([R sumStats.rBlank{u}],[],'off');
        isVis = signrank(R(:));

        if ~isAct || ~isVis
            continue
        end

        count = count+1;
        dist.rPref(count) = sumStats.rPref(u,setIdx);
        dist.dsi(count) = sumStats.dsi(u,setIdx);
        dist.ldr(count) = sumStats.ldr(u,setIdx);
        
    end

    if plt == 1
        figure; hold on
        
        subplot(2,3,1); hold on;
        cdfplot(dist.rPref)
        
        subplot(2,3,5); hold on;
        cdfplot(dist.dsi)
        
        subplot(2,3,6); hold on;
        cdfplot(dist.ldr)
    end
    
end