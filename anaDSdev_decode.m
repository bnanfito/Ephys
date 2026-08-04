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
nBoot = 100;
nFold = 5;
if strcmp(anaMode,'SU')
    nU = 34;
elseif strcmp(anaMode,'MU')
    nU = 300;
end
nPC = 30;

T2toT1 = @(in) rad2deg( mod(atan2(in(:,2),in(:,1)),2*pi) );
calcErr = @(in1,in2) rad2deg( angdiff(deg2rad( [in1,in2]' )) );

for ar = 1:nAR
for ag = 1:nAG

for b = 1:nBoot
    
    R = rMat{ar,ag}.rTrial_norm;
    R = R(:,randperm(size(R,2),nU));
    
    T1_dir = rMat{ar,ag}.cTrial;
    T1_dirShuff = T1_dir(randperm(length(T1_dir),length(T1_dir)));
    T1_ori = mod(T1_dir,180);
    T1_oriShuff = T1_ori(randperm(length(T1_ori),length(T1_ori)));

    T2_dir = [cos(deg2rad(T1_dir)) sin(deg2rad(T1_dir))]; %
    T2_dirShuff = [cos(deg2rad(T1_dirShuff)) sin(deg2rad(T1_dirShuff))];
    T2_ori = [cos(2*deg2rad(T1_ori)) sin(2*deg2rad(T1_ori))]; % multiply by 2 so that modulo 180 maps onto a unit circle instead of a semi-circle
    T2_oriShuff = [cos(deg2rad(T1_oriShuff)) sin(deg2rad(T1_oriShuff))];

    for p = 1:nPC
        pcs = 1:p;
        mae_dir = nan(1,nFold);
        mae_dirShuff = nan(1,nFold);
        mae_ori = nan(1,nFold);
        mae_oriShuff = nan(1,nFold);

        for f = 1:nFold
        
            cv = cvpartition(size(R,1),"KFold",nFold);
            trainIdx = training(cv,f);
            testIdx = test(cv,f);
            [cof,scr_train,~,~,~,mu] = pca(R(trainIdx,:));
            scr_test = ((R(testIdx,:)-mu)*cof(:,pcs));

            B_dir = scr_train(:,pcs)\T2_dir(trainIdx,:);
            T2pred_dir = scr_test*B_dir;
            T1pred_dir = T2toT1(T2pred_dir);
            err_dir = calcErr(T1_dir(testIdx),T1pred_dir);
            mae_dir(f) = mean(abs(err_dir));

            B_dirShuff = scr_train(:,pcs)\T2_dirShuff(trainIdx,:);
            T2pred_dirShuff = scr_test*B_dirShuff;
            T1pred_dirShuff = T2toT1(T2pred_dirShuff);
            err_dirShuff = calcErr(T1_dirShuff(testIdx),T1pred_dirShuff);
            mae_dirShuff(f) = mean(abs(err_dirShuff));
    
            B_ori = scr_train(:,pcs)\T2_ori(trainIdx,:);
            T2pred_ori = scr_test*B_ori;
            T1pred_ori = T2toT1(T2pred_ori);
            err_ori = calcErr(2*T1_ori(testIdx),T1pred_ori)/2;
            mae_ori(f) = mean(abs(err_ori));

            B_oriShuff = scr_train(:,pcs)\T2_oriShuff(trainIdx,:);
            T2pred_oriShuff = scr_test*B_oriShuff;
            T1pred_oriShuff = T2toT1(T2pred_oriShuff);
            err_oriShuff = calcErr(2*T1_oriShuff(testIdx),T1pred_oriShuff)/2;
            mae_oriShuff(f) = mean(abs(err_oriShuff));
        
%             if b==1 && p==1 && f==1
%                 figure
%                 plot(T2_dir(:,2),T2_dir(:,1),'ko');hold on
%                 plot([T2_dir(testIdx,2) T2pred_dir(:,2)]',[T2_dir(testIdx,1) T2pred_dir(:,1)]','r-')
%                 plot(T2pred_dir(:,2),T2pred_dir(:,1),'r*')
%                 axis square
%                 
%                 figure
%                 plot(T2_ori(:,2),T2_ori(:,1),'ko');hold on
%                 plot(T2pred_ori(:,2),T2pred_ori(:,1),'r*')
%                 plot([T2_ori(testIdx,2) T2pred_ori(:,2)]',[T2_ori(testIdx,1) T2pred_ori(:,1)]','r-')
%                 axis square
%             end

        end

        MAE_dir{ar,ag}(b,p) = mean(mae_dir); % compute mean mae across folds    
        MAE_dirShuff{ar,ag}(b,p) = mean(mae_dirShuff);
        MAE_ori{ar,ag}(b,p) = mean(mae_ori);
        MAE_oriShuff{ar,ag}(b,p) = mean(mae_oriShuff);

    end

end

end
end

%% plotting


figure; 

for ar = 1:nAR
for ag = 1:nAG
    subplot(nAR,nAG,ag+(nAG*(ar-1))); hold on
    
    yyaxis left
    y = mean(MAE_dir{ar,ag});
    err = confInt(MAE_dir{ar,ag});
    plot(y,'-o')
    patch([1:nPC fliplr(1:nPC)],[err(1,:) fliplr(err(2,:))],'k','EdgeColor','none','FaceAlpha',0.2)
    y = mean(MAE_dirShuff{ar,ag});
    err = confInt(MAE_dirShuff{ar,ag});
    plot(y)
    patch([1:nPC fliplr(1:nPC)],[err(1,:) fliplr(err(2,:))],'r','EdgeColor','none','FaceAlpha',0.2)
    ylim([40 120])
    xlim([0 nPC])
    
    yyaxis right
    y = mean(MAE_ori{ar,ag});
    err = confInt(MAE_ori{ar,ag});
    plot(y,'-o')
    patch([1:nPC fliplr(1:nPC)],[err(1,:) fliplr(err(2,:))],'k','EdgeColor','none','FaceAlpha',0.2)
    y = mean(MAE_oriShuff{ar,ag});
    err = confInt(MAE_oriShuff{ar,ag});
    plot(y)
    patch([1:nPC fliplr(1:nPC)],[err(1,:) fliplr(err(2,:))],'r','EdgeColor','none','FaceAlpha',0.2)
    ylim([20 60])
    xlim([0 nPC])
end
end