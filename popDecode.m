%popDecode

function [dis,guess,truth,confusion] = popDecode(r,c)

conds = unique(c);

nTrials = length(c);
nConds = length(unique(c));
nReps = nTrials/nConds;

% %scramble trials
% scrmbl = randperm(nTrials);
% [~,unscrmbl] = sort(scrmbl);
% r = r(scrmbl,:);
% c = c(scrmbl);

    for fold = 1:nReps
        testIdx = ismember(1:nTrials,fold:nReps:nTrials);
        trainIdx = ~testIdx;

        r_test = r(testIdx,:);
        c_test = c(testIdx);
        r_train = r(trainIdx,:);
        c_train = c(trainIdx);
        for i = 1:nConds
            f_train{fold}(i,:) = mean(r_train(c_train==conds(i),:),'omitnan');
        end

        for obs = 1:size(r_test,1)
            dis(obs,:,fold) = sum((r_test(obs,:)-f_train{fold}).^2,2);
            truth(obs,fold) = c_test(obs);
            minDisIdx = find(dis(obs,:,fold)==min(dis(obs,:,fold)));
            if length(minDisIdx)>1
                minDisIdx = minDisIdx(randi(length(minDisIdx)));
            end
            guess(obs,fold) = conds(minDisIdx);

        end

%         for i = 1:nConds
%             f_cur = f_train{fold}(i,:);
% 
%             dis(fold,i) = sum((r_test-f_cur).^2);
%         end
%         truth(fold) = c_test;
%         minC = conds(dis(fold,:)==min(dis(fold,:)));
%         if length(minC) == 1
%             guess(fold) = minC;
%         elseif length(minC) > 1
%             rndMin = randi(length(minC));
%             guess(fold) = minC(rndMin);
%         end
%         correct(fold) = guess(fold) == truth(fold);
    end
%     dis = dis(unscrmbl,:);
%     guess = guess(unscrmbl);
%     truth = truth(unscrmbl);
%     acc = sum(correct)/length(correct);

    for c = 1:nConds
        for g = 1:nConds
            confusion(g,c) = sum(guess(c,:)==conds(g));
        end
    end

%     figure;
%     subplot(2,2,1); hold on
%     imagesc(r)
%     axis tight
%     axis square
%     
%     subplot(2,2,2); hold on
%     imagesc(mean(dis,3))
%     cb = colorbar;
%     cb.Label.String = 'distance to condition mean';
%     axis tight
%     axis square
% 
%     subplot(2,2,3); hold on
%     imagesc()


end


















% figure;hold on
% for ag = 1:nAG
% 
%     r = rTrial{ag};
%     c = cTrial;
% %scramble trials
% scrmbl = randperm(nTrials);
% [~,unscrmbl] = sort(scrmbl);
% r = r(scrmbl,:);
% c = c(scrmbl);
%     subplot(2,nAG,ag)
%     imagesc(r(unscrmbl,:))
% 
%     for fold = 1:nTrials
%         testIdx = 1:nTrials == fold;
%         trainIdx = ~testIdx;
% 
%         r_test = r(testIdx,:);
%         c_test = c(testIdx);
%         r_train = r(trainIdx,:);
%         c_train = c(trainIdx);
%         for i = 1:nConds
%             f_train{fold,ag}(i,:) = mean(r_train(c_train==cMean(i),:),'omitnan');
%         end
% 
%         for i = 1:nConds
%             f_cur = f_train{fold,ag}(i,:);
% 
%             dis{ag}(fold,i) = sum((r_test-f_cur).^2);
%         end
%         truth{ag}(fold) = c_test;
%         guess{ag}(fold) = cMean(dis{ag}(fold,:)==min(dis{ag}(fold,:)));
%         correct{ag}(fold) = guess{ag}(fold) == truth{ag}(fold);
%     end
% 
%     acc(ag) = sum(correct{ag})/length(correct{ag});
%     subplot(2,nAG,ag+nAG)
%     imagesc(dis{ag}(unscrmbl,:))
%     axis tight
%     title(['age group ' num2str(ag) ': ' num2str(sum(correct{ag})) '/' num2str(length(correct{ag})) ' (' num2str(acc(ag)) '%)'])
% 
% end














% figure; hold on
% for ag = 1:nAG
% 
%     r = rTrial{ag};
%     c = cTrial;
% %scramble trials
% scrmbl = randperm(nTrials);
% [~,unscrmbl] = sort(scrmbl);
% r = r(scrmbl,:);
% c = c(scrmbl);
% 
%     for fold = 1:nReps
% 
%         testIdx = zeros(1,nTrials); testIdx(fold:nReps:nTrials) = 1; testIdx = testIdx == 1;
%         trainIdx = ~testIdx;
% 
% 
%         r_test = r(testIdx,:);
%         r_train = r(trainIdx,:);
%         c_train = c(trainIdx);
%         c_test = c(testIdx);
%         for i = 1:nConds
%             f_train{fold,ag}(i,:) = mean(r_train(c_train==cMean(i),:),'omitnan');
%         end
% 
%         nTest = size(r_test,1);
%         for t = 1:nTest
%             r_cur = r_test(t,:);
% 
%             for i = 1:nConds
%                 f_cur = f_train{fold,ag}(i,:);
% 
%                 dis{fold,ag}(t,i) = sum((r_cur-f_cur).^2);
%             end
%             guess{fold,ag}(t) = cMean(dis{fold,ag}(t,:)==min(dis{fold,ag}(t,:)));
%         end
%         correct = guess{fold,ag} == c_test;
%         acc(fold,ag) = sum(correct)/length(correct);
% 
%         subplot(nAG,nReps,fold+(nReps*(ag-1)))
%         imagesc(dis{fold,ag})
%     end
% 
% end
% 
% figure;hold on
% plot(1:nAG,acc,'k.')
% plot(1:nAG,mean(acc),'ro')