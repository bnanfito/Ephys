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

ar = 2;
ag = 1;

R = rMat{ar,ag}.rTrial_norm;
T1 = rMat{ar,ag}.cTrial;
T2 = [cos(deg2rad(T1)) sin(deg2rad(T1))];

nBoot = 100;
nK = 5;
nPC = 20;
for p = 1:nPC
    mae_boot = nan(nBoot,1);
    for b = 1:nBoot
        mae = nan(nK,1);
        for f = 1:nK
        
            cv = cvpartition(size(R,1),"KFold",nK);
            trainIdx = training(cv,f);
            testIdx = test(cv,f);
            
            [cof,scr,~,~,~,mu] = pca(R(trainIdx,:));
            
            pcs = p;
            B = scr(:,pcs)\T2(trainIdx,:);
            Tpred = ((R(testIdx,:)-mu)*cof(:,pcs))*B;
            Ttrain = ((R(trainIdx,:)-mu)*cof(:,pcs))*B;
            Tpred = rad2deg( mod(atan2(Tpred(:,2),Tpred(:,1)),2*pi) );
            Ttrain = rad2deg( mod(atan2(Ttrain(:,2),Ttrain(:,1)),2*pi) );
            
            err = rad2deg(angdiff(deg2rad(T1(testIdx)),deg2rad(Tpred)));
            mae(f) = mean(abs(err));
        
        end
        mae_boot(b) = mean(mae);
    end
    MAE(p) = mean(mae_boot);
    MAE_ci(:,p) = confInt(mae_boot);
end

% %plotting
% polaraxes
% 
% theta = repmat(deg2rad(T1(testIdx)),1,2)';
% rho = [zeros(size(theta,2),1) ones(size(theta,2),1)]';
% polarplot(theta,rho,'k-o');hold on
% 
% theta = repmat(deg2rad(T1(trainIdx)),1,2)';
% rho = [zeros(size(theta,2),1) ones(size(theta,2),1)]';
% polarplot(theta,rho,'b-o')
% 
% theta = [deg2rad(T1(trainIdx)) deg2rad(Ttrain)]';
% rho = ones(size(theta));
% polarplot(theta,rho,'g--*')
% 
% theta = [deg2rad(T1(testIdx)) deg2rad(Tpred)]';
% rho = ones(size(theta));
% polarplot(theta,rho,'r-*')

plot(MAE)
patch([1:nPC fliplr(1:nPC)],[MAE_ci(1,:) fliplr(MAE_ci(2,:))],'k','EdgeColor','none','FaceAlpha',0.2)