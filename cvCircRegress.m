
function [out] = cvCircRegress(r,t,varargin)

    cum = 0;
    maxPC = 14;

    % stimulus angle theta (t): direction and orientation values
    tDir = t;
    tDirShuff = tDir(randperm(length(tDir),length(tDir)));
    tOri = mod(tDir,180);
    tOriShuff = tOri(randperm(length(tOri),length(tOri)));

    % embed the stimulus angle theta in 2D space
    tD = [cos(deg2rad(tDir)) sin(deg2rad(tDir))];
    tDshuff = [cos(deg2rad(tDirShuff)) sin(deg2rad(tDirShuff))];
    tO = [cos(2*deg2rad(tOri)) sin(2*deg2rad(tOri))];
    tOshuff = [cos(2*deg2rad(tOriShuff)) sin(2*deg2rad(tOriShuff))];
    
    % reduce dimensionality
    [coeff,score] = pca(r);
    for nPC = 1:maxPC
        
        if cum == 1
            PCs = 1:nPC;
            yLb = 0;
        else
            PCs = nPC;
            yLb = 40;
        end
        nFold = 5;
        cv = cvpartition(length(tDir),'KFold',5);    
        errDir = zeros(1,nFold);
        errDirShuff = zeros(1,nFold);
        errOri = zeros(1,nFold);
        errOriShuff = zeros(1,nFold);
        for fold = 1:nFold
        
            trainIdx = training(cv,fold);
            testIdx  = test(cv,fold);
            cDir = score(trainIdx,PCs) \ tD(trainIdx,:); % learn mapping coeffs
            cDirShuff = score(trainIdx,PCs) \ tDshuff(trainIdx,:);
            cOri = score(trainIdx,PCs) \ tO(trainIdx,:);
            cOriShuff = score(trainIdx,PCs) \ tOshuff(trainIdx,:);
            tD_pred = score(testIdx,PCs)*cDir; % project PCs into 2D embedding
            tDir_pred = rad2deg(mod(atan2(tD_pred(:,2),tD_pred(:,1)),2*pi));
            tDshuff_pred = score(testIdx,PCs)*cDirShuff;
            tDirShuff_pred = rad2deg(mod(atan2(tDshuff_pred(:,2),tDshuff_pred(:,1)),2*pi));
            tO_pred = score(testIdx,PCs)*cOri;
            tOri_pred = rad2deg(mod(atan2(tO_pred(:,2),tO_pred(:,1)),2*pi)); %this is still in a modulo 2*pi circle, must be divided by two later
            tOshuff_pred = score(testIdx,PCs)*cOriShuff;
            tOriShuff_pred = rad2deg(mod(atan2(tOshuff_pred(:,2),tOshuff_pred(:,1)),2*pi));

            % diffDir = diff([tDir(testIdx) tDir_pred]');
            % diffOri = diff([2*tOri(testIdx) tOri_pred]')/2;
            diffDir = rad2deg(angdiff(deg2rad([tDir(testIdx) tDir_pred])'));
            diffDirShuff = rad2deg(angdiff(deg2rad([tDir(testIdx) tDirShuff_pred])'));
            diffOri = rad2deg(angdiff(deg2rad([2*tOri(testIdx) tOri_pred])'))/2;
            diffOriShuff = rad2deg(angdiff(deg2rad([2*tOri(testIdx) tOriShuff_pred])'))/2;

            errDir(fold) = mean(abs(diffDir));
            errDirShuff(fold) = mean(abs(diffDirShuff));
            errOri(fold) = mean(abs(diffOri));
            errOriShuff(fold) = mean(abs(diffOriShuff));
        
        end
        
        mseDir(nPC) = mean(errDir);
        mseDirShuff(nPC) = mean(errDirShuff);
        mseOri(nPC) = mean(errOri);
        mseOriShuff(nPC) = mean(errOriShuff);

    
    end

    out.errDir = mseDir;
    out.errDirShuff = mseDirShuff;
    out.errOri = mseOri;
    out.errOriShuff = mseOriShuff;


end