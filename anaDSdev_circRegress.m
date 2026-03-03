clear all
close all

%% Load data

anaMode = 'MU';

% dataFold = '/Volumes/Lab drive/Brandon/data/dataSets/DSdev';
% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/DSdev';
dataFold = '/Volumes/NielsenHome2/Brandon/data/dataSets/DSdev';
% dataFold = 'F:\Brandon\data\dataSets\DSdev';
% dataFold = 'Y:\Brandon\data\dataSets\DSdev';
load(fullfile(dataFold,['DSdev_' anaMode 'dataSet.mat']))

areas = {'V1','PSS'};
ageGroups = {[29 32],[33 36],[37 max(projTbl.age)]};
% ageGroups = {[28 32],[29 33],[30 34],[31 35],[32 36],[33 37],[34 38],[35 39],[36 40],[37 41],[38 42],[39 43],[40 44],[41 300]};
nAG = length(ageGroups);
nAR = length(areas);


figure;

maxPC = 10;
for ar = [2,1]
    for ag = 1:nAG

        for boot = 1:100
            disp(boot)
            r = rMat{ar,ag}.rTrial_norm;
            r = r(:,randperm(size(rMat{ar,ag}.rTrial_norm,2),300));
            tDir = rMat{ar,ag}.cTrial;
            tOri = mod(tDir,180);
            tD = [cos(deg2rad(tDir)) sin(deg2rad(tDir))];
            tO = [cos(2*deg2rad(tOri)) sin(2*deg2rad(tOri))];
            [coeff,score] = pca(r);
            for nPC = 1:maxPC
                
                nFold = 5;
                cv = cvpartition(length(tDir),'KFold',5);    
                errDir = zeros(1,nFold);
                errOri = zeros(1,nFold);
                for fold = 1:nFold
                
                    trainIdx = training(cv,fold);
                    testIdx  = test(cv,fold);
                    cDir = score(trainIdx,nPC) \ tD(trainIdx,:); % learn mapping coeffs
                    cOri = score(trainIdx,nPC) \ tO(trainIdx,:);
                    tD_pred = score(testIdx,nPC)*cDir; % project PCs into 2D embedding
                    tDir_pred = rad2deg(mod(atan2(tD_pred(:,2),tD_pred(:,1)),2*pi));
                    tO_pred = score(testIdx,nPC)*cOri;
                    tOri_pred = rad2deg(mod(atan2(tO_pred(:,2),tO_pred(:,1)),2*pi)); %this is still in a modulo 360 circle, must be divided by two later
                    
                    % diffDir = diff([tDir(testIdx) tDir_pred]');
                    % diffOri = diff([2*tOri(testIdx) tOri_pred]')/2;
                    diffDir = rad2deg(angdiff(deg2rad([tDir(testIdx) tDir_pred])'));
                    diffOri = rad2deg(angdiff(deg2rad([2*tOri(testIdx) tOri_pred])'))/2;
    
                    errDir(fold) = mean(abs(diffDir));
                    errOri(fold) = mean(abs(diffOri));
                
                end
                
                mseDir{ar,ag}(boot,nPC) = mean(errDir);
                mseOri{ar,ag}(boot,nPC) = mean(errOri);

                if nPC == 1 && ag==3 && ar==1 && boot==1
                    subplot(2,2,1);hold on
                    plot(tD(:,1),tD(:,2),'o')
                    plot(tD_pred(:,1),tD_pred(:,2),'r*')
                    plot(tD(testIdx,1),tD(testIdx,2),'g*')
                    plot([tD(testIdx,1) tD_pred(:,1)]',[tD(testIdx,2) tD_pred(:,2)]','k')
                    [x,y] = pol2cart(repmat(deg2rad(tDir_pred),1,2)',[zeros(length(tDir_pred),1) ones(length(tDir_pred),1)]');
                    plot(x,y,'r')
                    axis square
                    title('direction')
                    subplot(2,2,2);hold on
                    plot(tO(:,1),tO(:,2),'o')
                    plot(tO_pred(:,1),tO_pred(:,2),'r*')
                    plot(tO(testIdx,1),tO(testIdx,2),'g*')
                    plot([tO(testIdx,1) tO_pred(:,1)]',[tO(testIdx,2) tO_pred(:,2)]','k')
                    [x,y] = pol2cart(repmat(deg2rad(tOri_pred),1,2)',[zeros(length(tOri_pred),1) ones(length(tOri_pred),1)]');
                    plot(x,y,'r')
                    axis square
                    title('orientation')
                    colororder({'k','g'})
                end
            
            end
        end

    end
end

for ag = 1:nAG
    subplot(2,3,nAG+ag);hold on
    yyaxis left
    p(1) = plot(mean(mseDir{1,ag}),'k--');
    errY = confInt(mseDir{1,ag});
    patch([1:maxPC fliplr(1:maxPC)],[errY(1,:) fliplr(errY(2,:))],'k','EdgeColor','none','FaceAlpha',0.2)
    p(2) = plot(mean(mseDir{2,ag}),'k-');
    errY = confInt(mseDir{2,ag});
    patch([1:maxPC fliplr(1:maxPC)],[errY(1,:) fliplr(errY(2,:))],'k','EdgeColor','none','FaceAlpha',0.2)
    ylabel('direction error (mod 360)')
    ylim([0 180])
    yyaxis right
    p(3) = plot(mean(mseOri{1,ag}),'g--');
    errY = confInt(mseOri{1,ag});
    patch([1:maxPC fliplr(1:maxPC)],[errY(1,:) fliplr(errY(2,:))],'g','EdgeColor','none','FaceAlpha',0.2)
    p(4) = plot(mean(mseOri{2,ag}),'g-');
    errY = confInt(mseOri{2,ag});
    patch([1:maxPC fliplr(1:maxPC)],[errY(1,:) fliplr(errY(2,:))],'g','EdgeColor','none','FaceAlpha',0.2)
    ylabel('orientation error (mod 180)')
    ylim([0 90])
    if ag == 1
        legend(p,{'V1 dir','PMLS dir','V1 ori','PMLS ori'})
    end
end
