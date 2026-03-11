% anaRSA
clear all
% close all

%% generate tuning curves

x = 0:30:330;
a1 = 1;
mu1 = 180;
std1 = 10;
tc1 = normpdf(x,mu1,std1);
mu2 = mu1;
std2 = std1;
shift = find(x == mod(mu1+180,360))-find(x==180);
tc2 = normpdf(x,mu2,std2);
tc2 = circshift(tc2,shift);
tc2 = tc1+tc2;
tc = tc2;
tc = repmat(tc,100,1);
for i = 2:size(tc,1)
    tc(i,:) = circshift(tc(i,:),randi(size(tc,2)));
end
tc = tc';
tc = tc./max(tc);

%% rdm

diss = pdist(tc,'spearman');
rdm = squareform(diss);

%% dimensionality reduction

% score = mdscale(rdm,3);
[coeff,score] = pca(zscore(tc));

%% plot

figure

subplot(2,2,1);hold on
plot([x x(1)+360] ,[tc(:,1);tc(1,1)])
xlim([0 360])
xticks([0 90 180 270])

subplot(2,2,2);hold on
plot([x x(1)+360],([tc;tc(1,:)]))
xlim([0 360])
xticks([0 90 180 270])

subplot(2,2,3);hold on
imagesc(rdm)
axis tight

subplot(2,2,4);hold on
plot3([score(:,1);score(1,1)],[score(:,2);score(1,2)],[score(:,3);score(1,3)],'o-')
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
