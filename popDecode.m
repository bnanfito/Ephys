%popDecode

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

                dis(i) = sum((r_cur-f_cur).^2);
            end
            guess(fold,t) = cMean(dis==min(dis));
        end
        correct = guess(fold,:) == c_test;
        acc(fold,ag) = sum(correct)/length(correct);

        clear f_train dis guess
    end

end

figure;hold on
plot(1:nAG,acc,'k.')
plot(1:nAG,mean(acc),'ro')