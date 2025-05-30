function [cndInclude] = limitConds(sumStats)

    % if stimulus size is a variable being changes, only consider
    % conditions that were used to calculate sumStats (ie the ones that
    % were used to make the response matrix, R)
    u = 1; % input should have only one row (one unit)
    sizeId = contains(sumStats.paramKey{u},'size');
    cndSize = sumStats.cndKey{u}(:,sizeId);
    sizes = unique(cndSize);
    oriId = contains(sumStats.paramKey{u},'ori');
    cndOri = sumStats.cndKey{u}(:,oriId);
    conId = contains(sumStats.paramKey{u},'contrast');
    if sum(conId)>0
        R = sumStats.response{u}(:,:,sumStats.oriPref(u));
    else
        R = sumStats.response{u};
    end
    if sum(sizeId)>0 && length(sizes)>1
        if sum(conId)>0
            idx = cndSize>100 & cndOri==unique(sumStats.condition{u}(oriId,:,sumStats.oriPref(u)));
        else
            idx = cndSize>100;
        end
        ffR = sumStats(u,:).fr.bc(:,idx); %subset of FR mat for fullfield conditions
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

end