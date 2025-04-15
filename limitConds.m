function [cndInclude] = limitConds(sumStats)

    % if stimulus size is a variable being changes, only consider
    % conditions that were used to calculate sumStats (ie the ones that
    % were used to make the response matrix, R)
    u = 1; % input should have only one row (one unit)
    R = sumStats.response{u};
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

end