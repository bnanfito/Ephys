

function [R] = compute_rMat(data)

    nU = height(data);

    % %shuffle preferences
    % for u = 1:nU
    %     k = randi(size(data.response{u},2),1);
    %     data.response{u} = circshift(data.response{u},k,2);
    % end

    [~,oriPrefIdx] = sort(data.oriPref);
    data = data(oriPrefIdx,:); %sort units by their pref dir of motion
    
    % keepIdx = [];
    % for u = 1:nU
    %     cMean = sumStats.condition{u}(strcmp(sumStats.paramKey{u},'ori'),:);
    %     if length(cMean)==12
    %         keepIdx = [keepIdx u];
    %     end    
    % end
    % sumStats = sumStats(keepIdx,:);
    % nU = height(sumStats);
    
    for u = 1:nU
        cMean = data.condition{u}(strcmp(data.paramKey{u},'ori'),:);
        if length(cMean)==12
            rTrial(:,:,u) = data.response{u}(1:5,:);
        else
            cInt = 0:30:330;
            rTmp = data.response{u}(1:5,:);
            rInt = interp1([cMean cMean(1)+360],[rTmp rTmp(:,1)]',cInt);
            rTrial(:,:,u) = rInt';
            cMean = cInt;
        end
        rTrial(rTrial<0)=0;
        rMean(:,u) = mean(rTrial(:,:,u),'omitnan');
    
        cPref(1,u) = data.oriPref(u);
        cPref(2,u) = data.meanVec{u}.angDir;
    end
    nReps = size(rTrial,1);
    nConds = size(rTrial,2);
    rTrial = reshape(rTrial,nReps*nConds,nU);
    nTrial = size(rTrial,1);
    cTrial = repmat(cMean,nReps,1);cTrial = cTrial(:)';
    
    %standardize data
    rTrial_norm = rTrial./max(rTrial);
    rMean_norm = rMean./max(rMean);
    rTrial_z = zscore(rTrial);
    rMean_z = zscore(rMean);

    for u = 1:nU
        [tCent,rCent(:,u)] = alignDirTuning(cMean,rMean_norm(:,u)');
    end

    R.cPref = cPref;
    R.cMean = cMean';
    R.rMean = rMean;
    R.rMean_norm = rMean_norm;
    R.rMean_z = rMean_z;
    R.cTrial = cTrial';
    R.rTrial = rTrial;
    R.rTrial_norm = rTrial_norm;
    R.rTrial_z = rTrial_z;
    R.tCent = tCent;
    R.rCent = rCent;


end