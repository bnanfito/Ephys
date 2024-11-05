clear all
% close all

anaMode = 'SU';
proj = 'Train_V1Cool';
dataFold = 'D:\data\dataSets\training';
path = fullfile(dataFold,proj,anaMode,'ranksum & rPref above 2');
% path = fullfile(dataFold,proj,anaMode);
dataFileName = [proj '_' anaMode 'dataSet.mat'];
load(fullfile(path,dataFileName))

area = 'PSS';
exp = 'af';
ttl = [proj ': ' area exp ' training ' anaMode ' response'];
figFileName = [proj '_' anaMode '_' area exp '_popPCA.fig'];
if strcmp(area,'V1') && strcmp(exp,'bf')
    tbl = data.v1bf;
elseif strcmp(area,'V1') && strcmp(exp,'af')
    tbl = data.v1af;
elseif strcmp(area,'PSS') && strcmp(exp,'bf')
    tbl = data.pssbf;
elseif strcmp(area,'PSS') && strcmp(exp,'af')
    tbl = data.pssaf;
end
tbl = tbl(tbl.goodUnit,:);

% nU = height(tbl);
% for u = 1:nU
% 
%     rTemp = tbl.response{u};
%     cTemp = repmat(tbl.condition{u}(strcmp(tbl.paramKey{u},'ori'),:),size(rTemp,1),1);
% 
%     R(:,u) = rTemp(:);
%     C(:,u) = cTemp(:);
% end
% Rnorm = R./max(R);
% 
% r = Rnorm;
% c = mean(C,2); %collapse ori across units to be a 1d vectors (indexed by trial)
% 
% oris = unique(c);
% [coeff, score, latent, tsq, explained] = pca(r);
% 
% nClr = length(oris);
% clrs = hsv(nClr);
% figure; hold on
% 
% subplot(3,2,1);hold on
% imagesc(r)
% axis tight
% 
% subplot(3,2,2);hold on
% nP = size(score,1);
% for i = 1:nP
%     clr = clrs(oris == c(i),:);
%     plot3(score(i,1),score(i,2),score(i,3),'.','MarkerSize',20,'Color',clr)
% end
% % plot3(score(:,1),score(:,2),score(:,3),'k--')
% 
% subplot(3,2,3); hold on
% y = mean(r,2,'omitnan');
% plot(c,y)
% 
% subplot(3,2,4); hold on
% stem(cumsum(explained))
% xlabel('PC')
% ylabel('cum.sum explained variance')
% ylim([0 100])
% 
% subplot(3,1,3);hold on
% plot(c,score(:,1:4),'.','MarkerSize',20)



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
D = squareform(pdist(r));

for u = 1:nU
    [cAlign,rAlign(:,u),~] = alignDirTuning(C,r(:,u)');
end

figure; hold on
subplot(3,2,1);hold on
imagesc(r)
yticks(1:length(C))
yticklabels(num2str(C'))
xlabel('unit')
ylabel('direction (deg.)')
axis tight
title('norm. mean response matrix')

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

subplot(3,2,5);hold on
imagesc(D)
axis tight
xticks([1:size(D,2)])
xticklabels(num2str(C'))
yticks([1:size(D,1)])
yticklabels(num2str(C'))

subplot(3,2,6)
plot(C,score(:,1:4),'LineWidth',2)
% polarplot(deg2rad(C),score(:,1:4),'LineWidth',2)
% rlim([-1.5 1.5])
legend({'PC1','PC2','PC3','PC4'})

sgtitle(ttl)
saveas(gcf,fullfile(path,figFileName))

