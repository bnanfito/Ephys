%population vector

clear all
close all

proj = 'DSdev';
dataFold = '/Volumes/Lab drive/Brandon/data';

load(fullfile(dataFold,'dataSets',proj,'rMean_AG_pss.mat'))
rMean = r;
cMean = c;
clear r c

load(fullfile(dataFold,'dataSets',proj,'rTrial_AG_pss.mat'))
rTrial = r;
cTrial = c;
clear r c

nAG = length(rMean);
nTrials = length(cTrial);
nConds = length(cMean);
nReps = nTrials/nConds;
for ag = 1:nAG

    r = rTrial{ag};
    f = rMean{ag};

    for fold = 1:nReps

        testIdx = zeros(1,nTrials); testIdx(fold:nReps:nTrials) = 1; testIdx = testIdx == 1;
        trainIdx = ~testIdx;


        r_test = r(testIdx,:);
        r_train = r(trainIdx,:);
        c_train = cTrial(trainIdx);
        c_test = cTrial(testIdx);
        for i = 1:nConds
            f_train(i,:) = mean(r_train(c_train==cMean(i),:),'omitnan');
        end

        nTest = size(r_test,1);
        for t = 1:nTest
            r_cur = r_test(t,:);

            for i = 1:nConds
                f_cur = f_train(i,:);

%                 dis(i) = sum((r_cur-f_cur).^2);
            end
            guess(fold,t) = cMean(dis==min(dis));
        end
        correct = guess(fold,:) == c_test;
        acc(fold,ag) = sum(correct)/length(correct);

        clear f_train dis guess
    end

end



































% clear all
% close all
% 
% anaMode = 'MU';
% trainProj = 'Train_V1Cool';
% % dataFold = 'F:\Brandon\data';
% dataFold = '/Volumes/Lab drive/Brandon/data';
% dsFold = fullfile(dataFold,'dataSets','training',trainProj,anaMode);
% load(fullfile(dsFold,'ranksum & rPref above 2',[trainProj '_' anaMode 'dataSet.mat']))
% 
% stim = 180;
% for i = 1:4
% 
%     if i==1
%         I = 'v1bf';
%         tbl = data.v1bf;
%         clr = 'b';
%         linStyl = '--';
%     elseif i==2
%         I = 'v1af';
%         tbl = data.v1af;
%         clr = 'b';
%         linStyl = '-';
%     elseif i==3
%         I = 'pssbf';
%         tbl = data.pssbf;
%         clr = 'r';
%         linStyl = '--';
%     elseif i==4
%         I = 'pssaf';
%         tbl = data.pssaf;
%         clr = 'r';
%         linStyl = '-';
%     end
%     
%     nU = height(tbl);
%     for u = 1:nU
%     
%         angDir(u) = tbl.meanVec{u}.angDir;
%         ldr(u) = tbl.meanVec{u}.ldr;
%         oriInd = find(strcmp(tbl.paramKey{u},'ori'));
%         tune(u,:) = mean(tbl.response{u},'omitnan');
%         ori(u,:) = tbl.condition{u}(oriInd,:);
%         
%     
%     %     polarplot(repmat(deg2rad(angDir),1,2),[0 ldr],'k');hold on
%     
%     end
%     
%     w = tune(ori==stim);
%     pV = sum(w.*exp(sqrt(-1)*mod(deg2rad(angDir(:)),2*pi)));
%     
%     
%     PV(i) = rad2deg(mod( angle(pV) ,2*pi));
% 
% 
% end
