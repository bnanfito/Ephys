%popDecode

clear all
close all

proj = 'DSdev';
dataFold = '/Volumes/Lab drive/Brandon/data';

load(fullfile(dataFold,'dataSets',proj,'rMean_AG.mat'))
rMean = r;
cMean = c;
clear r c

load(fullfile(dataFold,'dataSets',proj,'rTrial_AG.mat'))
rTrial = r;
cTrial = c;
clear r c

nAG = length(rMean);
nTrials = length(cTrial);
nConds = length(cMean);
nReps = nTrials/nConds;
for ag = 2

    r = rTrial{ag};
    f = rMean{ag};

    for fold = 1:5

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
            guess(t) = cMean(dis==min(dis));
        end
        correct = guess == c_test;
        acc(fold) = sum(correct)/length(correct);

        clear f_train dis guess
    end

end