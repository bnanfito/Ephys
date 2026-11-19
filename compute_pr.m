% compute_pr.m
% compute participation ratio for estimation of linear dimensionality

function [D] = compute_pr(R)

cv = 1;
if cv==1
    
    for b = 1:100
        idx = randi(2,size(R,1),1);
        A = R(idx==1,:);
        B = R(idx==2,:);
        [P] = pca(B);
        C = cov(A);
        for i = 1:size(P,2)
            V(b,i) = (P(:,i)'*C*P(:,i))/trace(C);
        end

    end
    L = mean(V);
    D = (sum(L)^2)/sum(L.^2);

elseif cv==0

    [~,~,L] = pca(R); 
    D = (sum(L)^2)/sum(L.^2);

%     C = cov(R);
%     D = (trace(C)^2)/trace(C^2);

end

% nBoot = 100;
% nC = size(R,2);
% D = nan(nBoot,nC);
% for n = 1:5:nC
%     for b = 1:nBoot
%         idx = randperm(nC,n);
%         Rn = R(:,idx);
%         
%         C = cov(Rn);
%         D(b,n) = (trace(C)^2) / trace(C^2);
%     end
% end
% 
% figure(1);hold on
% stem(mean(D))
% ci = confInt(D);
% patch([1:nC fliplr(1:nC)],...
%     [ci(1,:) fliplr(ci(2,:))],...
%     'r','FaceAlpha',0.2)

end