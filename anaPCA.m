clear all
close all

anaMode = 'SU';
proj = 'Train_V1Cool';
path = ['/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/training/' proj '/' anaMode '/ranksum & rPref above 2/'];
dataFileName = [proj '_' anaMode 'dataSet.mat'];
load([path dataFileName])

ttl = [proj ': PSS after training ' anaMode ' response'];
figFileName = [proj '_' anaMode 'pssaf_popPCA.fig'];
tbl = data.pssaf;
tbl = tbl(tbl.goodUnit,:);

% [~,lateSortIdx] = sort(tbl.latency);
% tBin = 0.01;
% for u = 1:height(tbl)
%     x = -1:tBin:2;
%     spks = tbl.spkTimes{u}(1,:);
%     g = normpdf(repmat(x,length(spks),1)',spks,repmat(0.05,1,length(spks)));
%     sdf(u,:) = mean(g,2,'omitnan');
%     sdf(u,:) = sdf(u,:)/max(sdf(u,:));
% end
% 
% figure; hold on
% imagesc(sdf(lateSortIdx,:))
% ylabel('unit')
% xlabel(['time bin (' num2str(tBin) ' sec)'])
% title('SDF')
% axis tight


R = cat(3,tbl.response{:});
R(R<0) = 0;
nU = size(R,3);
for u = 1:nU
%     Z(u,:) = ()
end
Rmean = squeeze(mean(R,1,'omitnan'));
Rnorm = Rmean./max(Rmean);
C = tbl.condition{1};

r = Rnorm;
[coeff, score, latent, tsq, explained] = pca(r);

for u = 1:nU
    [cAlign,rAlign(:,u),~] = alignDirTuning(C,r(:,u)');
end


figure; hold on
subplot(3,2,1);hold on
imagesc(r')
xticks(1:length(C))
xticklabels(mat2cell(C,1,ones(1,length(C))))
ylabel('unit')
xlabel('direction (deg.)')
axis tight
colormap(gca)
title('norm. response matrix')


subplot(3,2,3);hold on
plot(cAlign,rAlign)
plot(cAlign,mean(rAlign,2),'k','LineWidth',2)
sem = std(rAlign,[],2)/sqrt(size(rAlign,2));
plot(repmat(cAlign,2,1),mean(rAlign,2)'+([1;-1]*sem'),'k','LineWidth',2)
xlabel('direction (deg. from pref.)')
ylabel('response (normalized)')
title('tuning curve')


subplot(3,2,2);hold on
nP = size(score,1);
clrs = hsv(nP);
for i = 1:nP
    pt(i) = plot3(score(i,1),score(i,2),score(i,3),'.','Color',clrs(i,:),'MarkerSize',20);
end
plot3(score(:,1),score(:,2),score(:,3),'k--')
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
legend(pt,num2str(C'))
title('PCA')

subplot(3,2,4);hold on
stem(cumsum(explained))
xlabel('PC')
ylabel('cum.sum explained variance')
ylim([0 100])

subplot(3,1,3)
plot(C,score(:,1:4),'LineWidth',2)
% polarplot(deg2rad(C),score(:,1:4),'LineWidth',2)
% rlim([-1.5 1.5])
legend({'PC1','PC2','PC3','PC4'})

sgtitle(ttl)
% saveas(gcf,[path figFileName])

