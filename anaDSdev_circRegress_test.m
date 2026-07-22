clear all
close all

%% Load data

anaMode = 'SU';

% dataFold = '/Volumes/Lab drive/Brandon/data/dataSets/DSdev';
% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/DSdev';
% dataFold = '/Volumes/NielsenHome2/Brandon/data/dataSets/DSdev';
% dataFold = 'F:\Brandon\data\dataSets\DSdev';
dataFold = 'Y:\Brandon\data\dataSets\DSdev';
load(fullfile(dataFold,['DSdev_' anaMode 'dataSet.mat']))

areas = {'V1','PSS'};
ageGroups = {[29 32],[33 36],[37 max(projTbl.age)]};
% ageGroups = {[28 32],[29 33],[30 34],[31 35],[32 36],[33 37],[34 38],[35 39],[36 40],[37 41],[38 42],[39 43],[40 44],[41 300]};
nAG = length(ageGroups);
nAR = length(areas);


figure;
cum = 0;
exAR = 1;
exAG = 2;
maxPC = 10;
exPC = 1;
if strcmp(anaMode,'SU')
    nU = 34;
elseif strcmp(anaMode,'MU')
    nU = 300;
end
for ar = 1:nAR
    for ag = 1:nAG

        for boot = 1:100
            disp(boot)
            r = rMat{ar,ag}.rTrial_norm;
            r = r(:,randperm(size(r,2),nU));
            tDir = rMat{ar,ag}.cTrial;
            tDirShuff = tDir(randperm(length(tDir),length(tDir)));
            tOri = mod(tDir,180);
            tOriShuff = tOri(randperm(length(tOri),length(tOri)));
%             tD = [cos(deg2rad(tDir)) sin(deg2rad(tDir))];
%             tDshuff = [cos(deg2rad(tDirShuff)) sin(deg2rad(tDirShuff))];
%             tO = [cos(2*deg2rad(tOri)) sin(2*deg2rad(tOri))];
%             tOshuff = [cos(2*deg2rad(tOriShuff)) sin(2*deg2rad(tOriShuff))];
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
                    cDir = score(trainIdx,PCs) \ tDir(trainIdx,:); % learn mapping coeffs
                    cDirShuff = score(trainIdx,PCs) \ tDirShuff(trainIdx,:);
                    cOri = score(trainIdx,PCs) \ tOri(trainIdx,:);
                    cOriShuff = score(trainIdx,PCs) \ tOriShuff(trainIdx,:);
                    tDir_pred = score(testIdx,PCs)*cDir; % project PCs into 2D embedding
%                     tDir_pred = rad2deg(mod(atan2(tD_pred(:,2),tD_pred(:,1)),2*pi));
                    tDirShuff_pred = score(testIdx,PCs)*cDirShuff;
%                     tDirShuff_pred = rad2deg(mod(atan2(tDshuff_pred(:,2),tDshuff_pred(:,1)),2*pi));
                    tOri_pred = score(testIdx,PCs)*cOri;
%                     tOri_pred = rad2deg(mod(atan2(tO_pred(:,2),tO_pred(:,1)),2*pi)); %this is still in a modulo 2*pi circle, must be divided by two later
                    tOriShuff_pred = score(testIdx,PCs)*cOriShuff;
%                     tOriShuff_pred = rad2deg(mod(atan2(tOshuff_pred(:,2),tOshuff_pred(:,1)),2*pi));

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
                
                mseDir{ar,ag}(boot,nPC) = mean(errDir);
                mseDirShuff{ar,ag}(boot,nPC) = mean(errDirShuff);
                mseOri{ar,ag}(boot,nPC) = mean(errOri);
                mseOriShuff{ar,ag}(boot,nPC) = mean(errOriShuff);

                if nPC == exPC && ag==exAG && ar==exAR && boot==1
                    subplot(2,2,1,polaraxes)
                    polarplot(deg2rad(tDir),1,'o');hold on
                    polarplot( repmat(deg2rad(tDir(testIdx)),1,2)',[zeros(sum(testIdx),1) ones(sum(testIdx),1)]','k--','LineWidth',1)
                    polarplot(deg2rad(tDir_pred),ones(size(tDir_pred)),'r*')
%                     polarplot([tDir(testIdx,1) tDir_pred(:,1)]',[tD(testIdx,2) tD_pred(:,2)]','k')
%                     [x,y] = pol2cart(repmat(deg2rad(tDir_pred),1,2)',[zeros(length(tDir_pred),1) ones(length(tDir_pred),1)]');
                    polarplot(repmat(deg2rad(tDir_pred),1,2)',[zeros(length(tDir_pred),1) ones(length(tDir_pred),1)]','r','LineWidth',1)
%                     plot([x(2,:);tD(testIdx,1)'], [y(2,:);tD(testIdx,2)'],'b','LineWidth',1)
%                     axis square
%                     xlim([-1.2 1.2]);ylim([-1.2 1.2])
%                     box on
                    title('direction')
                    subplot(2,2,2,polaraxes)
                    polarplot(deg2rad(tOri),ones(size(tOri)),'o');hold on
                    polarplot( repmat(deg2rad(tOri(testIdx)),1,2)',[zeros(sum(testIdx),1) ones(sum(testIdx),1)]','k--','LineWidth',1)
                    polarplot(deg2rad(tOri_pred),ones(size(tOri_pred)),'r*')
%                     plot([tO(testIdx,1) tO_pred(:,1)]',[tO(testIdx,2) tO_pred(:,2)]','k')
%                     [x,y] = pol2cart(repmat(deg2rad(tOri_pred),1,2)',[zeros(length(tOri_pred),1) ones(length(tOri_pred),1)]');
                    polarplot(repmat(deg2rad(tOri_pred),1,2)',[zeros(length(tOri_pred),1) ones(length(tOri_pred),1)]','r','LineWidth',1)
%                     plot([x(2,:);tO(testIdx,1)'], [y(2,:);tO(testIdx,2)'],'g','LineWidth',1)
%                     axis square
%                     xlim([-1.2 1.2]);ylim([-1.2 1.2])
%                     box on
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
    p(1) = plot(mean(mseDir{1,ag}),'ko--','LineWidth',1);
    errY = confInt(mseDir{1,ag});
    patch([1:maxPC fliplr(1:maxPC)],[errY(1,:) fliplr(errY(2,:))],'k','EdgeColor','none','FaceAlpha',0.2)
    p(2) = plot(mean(mseDir{2,ag}),'ko-','LineWidth',1);
    errY = confInt(mseDir{2,ag});
    patch([1:maxPC fliplr(1:maxPC)],[errY(1,:) fliplr(errY(2,:))],'k','EdgeColor','none','FaceAlpha',0.2)
    plot(mean(mseDirShuff{1,ag}),'r--')
    errY = confInt(mseDirShuff{1,ag});
    patch([1:maxPC fliplr(1:maxPC)],[errY(1,:) fliplr(errY(2,:))],'r','EdgeColor','none','FaceAlpha',0.2)
    plot(mean(mseDirShuff{2,ag}),'r-')
    errY = confInt(mseDirShuff{2,ag});
    patch([1:maxPC fliplr(1:maxPC)],[errY(1,:) fliplr(errY(2,:))],'r','EdgeColor','none','FaceAlpha',0.2)
    if ag == 1
        ylabel('direction error (mod 360)')
    end
    ylim([yLb 120])
    yyaxis right
    p(3) = plot(mean(mseOri{1,ag}),'go--','LineWidth',1);
    errY = confInt(mseOri{1,ag});
    patch([1:maxPC fliplr(1:maxPC)],[errY(1,:) fliplr(errY(2,:))],'g','EdgeColor','none','FaceAlpha',0.2)
    p(4) = plot(mean(mseOri{2,ag}),'go-','LineWidth',1);
    errY = confInt(mseOri{2,ag});
    patch([1:maxPC fliplr(1:maxPC)],[errY(1,:) fliplr(errY(2,:))],'g','EdgeColor','none','FaceAlpha',0.2)
    p(5) = plot(mean(mseOriShuff{1,ag}),'r--');
    errY = confInt(mseOriShuff{1,ag});
    patch([1:maxPC fliplr(1:maxPC)],[errY(1,:) fliplr(errY(2,:))],'r','EdgeColor','none','FaceAlpha',0.2)
    plot(mean(mseOriShuff{2,ag}),'r-')
    errY = confInt(mseOriShuff{2,ag});
    patch([1:maxPC fliplr(1:maxPC)],[errY(1,:) fliplr(errY(2,:))],'r','EdgeColor','none','FaceAlpha',0.2)
    ylim([yLb/2 60])
    if ag == 3
        ylabel('orientation error (mod 180)')
    end
    if ag == 3
        legend(p,{'V1 dir','PMLS dir','V1 ori','PMLS ori','chance'})
    elseif ag == 2
        xlabel('number of PCs')
    end
    title(['AG' num2str(ag)])
    box on
end
sgtitle([areas{exAR} ' AG' num2str(exAG) ' ' anaMode ' circular regression example fold (PC#1-' num2str(exPC) ')'])


figure;
for ar = 1:nAR
    if strcmp(areas{ar},'V1')
        clrDir = 'b';
        clrOri = 'c';
    elseif strcmp(areas{ar},'PSS')
        clrDir = 'r';
        clrOri = 'm';
    end
    colororder = {clrDir clrOri};
    subplot(1,2,ar);hold on
    box on
    for ag = 1:nAG
        if ag == 1
            linStyl = ':';
            mrkr = 'o';
        elseif ag == 2
            linStyl = '--';
            mrkr = '^';
        elseif ag == 3
            linStyl = '-';
            mrkr = 'square';
        end
        yyaxis left
        plot(mean(mseDir{ar,ag}),[clrDir linStyl mrkr],'LineWidth',1);
        errY = confInt(mseDir{ar,ag});
        patch([1:maxPC fliplr(1:maxPC)],[errY(1,:) fliplr(errY(2,:))],clrDir,'EdgeColor','none','FaceAlpha',0.2)
        ylim([yLb 120])
        yyaxis right
        plot(mean(mseOri{ar,ag}),[clrOri linStyl mrkr],'LineWidth',1);
        errY = confInt(mseOri{ar,ag});
        patch([1:maxPC fliplr(1:maxPC)],[errY(1,:) fliplr(errY(2,:))],clrOri,'EdgeColor','none','FaceAlpha',0.2)
        ylim([yLb/2 60])
    end



end