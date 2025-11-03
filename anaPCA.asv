% close all
% clear

function [distDat] = anaPCA(sumStats)

plt = 0;
tAve = 0;

[~,oriPrefIdx] = sort(sumStats.oriPref);
sumStats = sumStats(oriPrefIdx,:); %sort units by their pref dir of motion

nU = height(sumStats);
for u = 1:nU
    cMean = sumStats.condition{u}(strcmp(sumStats.paramKey{1},'ori'),:);
    if length(cMean)==12
        rTrial(:,:,u) = sumStats.response{u}(1:5,:);
    else
        cInt = 0:30:330;
        rTmp = sumStats.response{u}(1:5,:);
        rInt = interp1([cMean cMean(1)+360],[rTmp rTmp(:,1)]',cInt);
        rTrial(:,:,u) = rInt';
        cMean = cInt;
    end
    rTrial(rTrial<0)=0;
    rMean(:,u) = mean(rTrial(:,:,u),'omitnan');

    cPref(1,u) = sumStats.oriPref(u);
    cPref(2,u) = sumStats.meanVec{u}.angDir;
end
nReps = size(rTrial,1);
nConds = size(rTrial,2);
rTrial = reshape(rTrial,nReps*nConds,nU);
nTrial = size(rTrial,1);
cTrial = repmat(cMean,nReps,1);cTrial = cTrial(:)';

%standardize data
rTrial_norm = rTrial./max(rTrial);
rMean_norm = rMean./max(rMean);
rTrial_z = zscore(rTrial);
rMean_z = zscore(rMean);


if tAve == 1
    x = rMean_norm;
    y = cMean;
else
    x = rTrial_norm;
    y = cTrial;
end
np = size(x,1);
if tAve == 1
    oriClrs = repmat(hsv(np/2),2,1);
    dirClrs = hsv(np);
    clrs = dirClrs;
else
    clrsTmp = hsv(nConds);
    clrs = [];
    for i = 1:nConds
        clrs = vertcat(clrs,repmat(clrsTmp(i,:),nReps,1));
    end
end

% PCA
[coeff,score,latent,tsquare,explained] = pca(x);

% RDM
diss = pdist(x,'spearman');
rdm = squareform(diss);
% rdm = 1-corr(x','type','Spearman');
% rdm = dist(x');
% rdm = squareform(pdist(x,'squaredeuclidean'));

for u = 1:nU
    [tCent,rCent(:,u)] = alignDirTuning(y,x(:,u)');
end

%% distance function

for i = 1:size(rdm,1)
    shift = size(rdm,1)-(i-1);
    rdmShift(i,:) = circshift(rdm(i,:),shift,2);
end
distF = mean(rdmShift);

%% Null Dist

nIter = 10000;
for iter = 1:nIter

    randIdx = randperm(size(x,1));
    shuff = x(randIdx,:);

    shuff = shuff./max(shuff);

    [coeffShuff,scoreShuff,latentShuff,tsquareShuff,explainedShuff] = pca(shuff);
    dissNull(:,iter) = pdist(shuff,'spearman');
    rdmNull(:,:,iter) = squareform(dissNull(:,iter));
%     rdmNull(:,:,iter) = 1-corr(shuff','type','Spearman');
%     rdmNull(:,:,nr) = dist(shuff');

    for i = 1:size(rdmNull,1)
        shift = size(rdmNull,1)-(i-1);
        rdmNullShift(i,:,iter) = circshift(rdmNull(i,:,iter),shift,2);
    end
    distNull(iter,:) = mean(rdmNullShift(:,:,iter),1,'omitnan');

end

distF_z = (distF-mean(distNull,'omitnan'))./std(distNull,'omitnan');

%% Output structure

distDat.cPref = cPref;
distDat.cMean = cMean';
distDat.rMean = rMean;
distDat.rMean_norm = rMean_norm;
distDat.rMean_z = rMean_z;
distDat.cTrial = cTrial';
distDat.rTrial = rTrial;
distDat.rTrial_norm = rTrial_norm;
distDat.rTrial_z = rTrial_z;
distDat.tCent = tCent;
distDat.rCent = rCent;
distDat.score = score;
distDat.coeff = coeff;
distDat.diss = diss;
distDat.rdm = rdm;
distDat.rdmShift = rdmShift;
distDat.distF = distF;
distDat.dissNull = dissNull;
distDat.rdmNull = rdmNull;
distDat.distNull = distNull;

%% Plot

if plt == 1



figure;hold on;
imagesc(x);
axis tight
if tAve == 1
    yticks(1:size(x,1))
    yticklabels(num2str(y'))
else
    yticks(3:nReps:size(x,1))
    yticklabels(num2str(y(1:nReps:size(x,1))'))
end
colorbar
xlabel('unit')
ylabel('direction of motion')



figure;hold on;
plot(tCent,rCent) 
plot(tCent,mean(rCent,2,'omitnan'),'k','LineWidth',2)
plot(repmat(tCent,2,1),mean(rCent,2,'omitnan')'+([1;-1].*(std(rCent,[],2,'omitnan')/sqrt(size(rCent,2)))'),'k','LineWidth',2)
xlabel('deg. relative to preferred')
ylabel('mean (+/-sem) normalized response')



figure;hold on;
imagesc(rdm);
axis tight
yticks(1:length(y))
yticklabels(num2str(y'))
xticks(1:length(y))
xticklabels(num2str(y'))
colorbar
xBox = [0 1 1 0]; xOrthOff = [(0.5:11.5) (0.5:3.5)]; xBox = repmat(xBox',1,length(xOrthOff))+xOrthOff;
yBox = [0 0 1 1]; yOrthOff = [(4.5:15.5) (12.5:15.5)]; yBox = repmat(yBox',1,length(yOrthOff))+yOrthOff;
% patch(yBox,xBox,'r','FaceColor','none','EdgeColor','r','LineWidth',2)
xlabel('direction of motion')
ylabel('direction of motion')
title('RDM')



figure;hold on;
if tAve == 1
    plot3([scoreShuff(:,1);scoreShuff(1,1)],[scoreShuff(:,2);scoreShuff(1,2)],[scoreShuff(:,3);scoreShuff(1,3)],'--','Color',[0.8 0.8 0.8])
end
for i = 1:np
    clr = clrs(i,:);
    pt(i) = plot3(score(i,1),score(i,2),score(i,3),'.','Color',clr,'MarkerSize',20);
end
if tAve == 1
    plot3([score(:,1);score(1,1)],[score(:,2);score(1,2)],[score(:,3);score(1,3)],'k--','LineWidth',2)
    legend(pt,num2str(y'))
end
xlabel('PC1');ylabel('PC2');zlabel('PC3')



figure;hold on
stem(cumsum(explained))
xlabel('PC')
ylabel('cum.sum explained variance')
ylim([0 100])



figure
if tAve == 1
    plot(y,score(:,1:4),'LineWidth',2)
else
    plot(y,score(:,1:4),'o')
end
% polarplot(deg2rad(cMean),score(:,1:4),'LineWidth',2)
% rlim([-1.5 1.5])
legend({'PC1','PC2','PC3','PC4'})



figure;hold on
dispIdx = y<=180;
angDisp = y(dispIdx);
distF = distF(dispIdx);
distNull = distNull(:,dispIdx);
% plot(angDisp,distF_z,'k-o','LineWidth',2)
plot(angDisp,distF,'k-o','LineWidth',2)
sem = std(rdmShift)/sqrt(size(rdmShift,1)); sem = sem(dispIdx);
v = var(rdmShift); v = v(dispIdx);
patch([angDisp fliplr(angDisp)],[distF-v fliplr(distF+v)],'k','EdgeColor','none','FaceAlpha',0.2)

plot(angDisp,mean(distNull),'r-o')
sem = std(distNull)/sqrt(size(distNull,1));
v = var(distNull);
sig2 = std(distNull,'omitnan')*2;
for i = 1:size(distNull,2)
    P99(i) = prctile(distNull(:,i),99,'Method','exact');
    P95(i) = prctile(distNull(:,i),95,'Method','exact');
    P05(i) = prctile(distNull(:,i),5,'Method','exact');
    P01(i) = prctile(distNull(:,i),1,'Method','exact');
end
patch([angDisp fliplr(angDisp)],[mean(distNull)-v fliplr(mean(distNull)+v)],'r','EdgeColor','none','FaceAlpha',0.2)
patch([angDisp fliplr(angDisp)],[P05 fliplr(P95)],'r','FaceColor','none','EdgeColor','r','LineStyle','--')
patch([angDisp fliplr(angDisp)],[P01 fliplr(P99)],'r','FaceColor','none','EdgeColor','r','LineStyle',':')

% plot(angDisp,mean(distNull2),'b--*')
% sem2 = std(distNull2)/sqrt(size(distNull2,1));
% v2 = var(distNull2);
% sig2_2 = std(distNull2,'omitnan')*2;
% for i = 1:size(distNull2,2)
%     h = cdfplot(distNull2(:,i));
%     P95_2(i) = h.XData(find(h.YData>=0.95,1));
%     P05_2(i) = h.XData(find(h.YData>=0.05,1));
%     delete(h)
% end
% patch([angDisp fliplr(angDisp)],[mean(distNull2)-v2 fliplr(mean(distNull2)+v2)],'b','EdgeColor','none','FaceAlpha',0.2)
% patch([angDisp fliplr(angDisp)],[P05_2 fliplr(P95_2)],'b','FaceColor','none','EdgeColor','b','LineStyle','--')

sigHiX = angDisp(distF>P95);
sigHiY = distF(distF>P95);
sigLoX = angDisp(distF<P05);
sigLoY = distF(distF<P05);
text(sigHiX,sigHiY+(sigHiY*0.1),'*')
text(sigLoX,sigLoY-(sigLoY*0.1),'*')

xticks([0 90 180])
xticklabels({'0','90','180'})
xlabel('angular disparity (+/- deg)')



end


end






















% clear all
% close all
% 
% anaMode = 'MU';
% proj = 'Train_MonoBiDir';
% dataFold = 'Y:\Brandon\data\dataSets\training';
% % path = fullfile(dataFold,proj,anaMode,'ranksum & rPref above 2');
% path = fullfile(dataFold,proj,anaMode);
% dataFileName = [proj '_' anaMode 'dataSet.mat'];
% load(fullfile(path,dataFileName))
% 
% area = 'PSS';
% exp = 'af';
% ttl = [proj ': ' area exp ' training ' anaMode ' response'];
% figFileName = [proj '_' anaMode '_' area exp '_popPCA.fig'];
% if strcmp(area,'V1') && strcmp(exp,'bf')
%     tbl = data.v1bf;
% elseif strcmp(area,'V1') && strcmp(exp,'af')
%     tbl = data.v1af;
% elseif strcmp(area,'PSS') && strcmp(exp,'bf')
%     tbl = data.pssbf;
% elseif strcmp(area,'PSS') && strcmp(exp,'af')
%     tbl = data.pssaf;
% end
% tbl = tbl(tbl.goodUnit,:);
% [~,sortIdx] = sort(tbl.oriPref);
% tbl = tbl(sortIdx,:);
% nU = height(tbl);
% 
% 
% 
% % nU = height(tbl);
% % for u = 1:nU
% % 
% %     rTemp = tbl.response{u};
% %     cTemp = repmat(tbl.condition{u}(strcmp(tbl.paramKey{u},'ori'),:),size(rTemp,1),1);
% % 
% %     R(:,u) = rTemp(:);
% %     C(:,u) = cTemp(:);
% % end
% % Rnorm = R./max(R);
% % 
% % r = Rnorm;
% % c = mean(C,2); %collapse ori across units to be a 1d vectors (indexed by trial)
% % 
% % oris = unique(c);
% % [coeff, score, latent, tsq, explained] = pca(r);
% % D = squareform(pdist(r));
% % 
% % nClr = length(oris);
% % clrs = hsv(nClr);
% % figure; hold on
% % 
% % subplot(3,2,1);hold on
% % imagesc(r)
% % axis tight
% % colorbar
% % ylabel('trial')
% % xlabel('unit')
% % 
% % subplot(3,2,2);hold on
% % nP = size(score,1);
% % for i = 1:nP
% %     clr = clrs(oris == c(i),:);
% %     plot3(score(i,1),score(i,2),score(i,3),'.','MarkerSize',20,'Color',clr)
% % end
% % xlabel('PC1')
% % ylabel('PC2')
% % zlabel('PC3')
% % 
% % subplot(3,2,3); hold on
% % y = mean(r,2,'omitnan');
% % plot(c,y)
% % 
% % subplot(3,2,4); hold on
% % stem(cumsum(explained))
% % xlabel('PC')
% % ylabel('cum.sum explained variance')
% % ylim([0 100])
% % 
% % subplot(3,2,5);hold on
% % imagesc(D)
% % axis tight
% % 
% % subplot(3,2,6);hold on
% % plot(c,score(:,1:4),'.','MarkerSize',20)
% % 
% % sgtitle(ttl)
% % saveas(gcf,fullfile(path,['trial_' figFileName]))
% 
% 
% 
% % R = cat(3,tbl.response{:});
% % R(R<0) = 0;
% % Rmean = squeeze(mean(R,1,'omitnan'));
% for u = 1:nU
%     R = tbl.response{u};
%     R(R<0) = 0;
%     Rmean(:,u) = mean(R,'omitnan');
% end
% Rnorm = Rmean./max(Rmean);
% C = tbl.condition{1}(strcmp(tbl.paramKey{1},'ori'),:);
% 
% r = Rnorm;
% [coeff, score, latent, tsq, explained] = pca(r);
% D = squareform(pdist(r));
% 
% for u = 1:nU
%     [cAlign,rAlign(:,u),~] = alignDirTuning(C,r(:,u)');
% end
% 
% figure; hold on
% subplot(3,2,1);hold on
% imagesc(r)
% yticks(1:length(C))
% yticklabels(num2str(C'))
% xlabel('unit')
% ylabel('direction (deg.)')
% axis tight
% title('norm. mean response matrix')
% colorbar
% 
% subplot(3,2,3);hold on
% plot(cAlign,rAlign)
% plot(cAlign,mean(rAlign,2),'k','LineWidth',2)
% sem = std(rAlign,[],2)/sqrt(size(rAlign,2));
% plot(repmat(cAlign,2,1),mean(rAlign,2)'+([1;-1]*sem'),'k','LineWidth',2)
% xlabel('direction (deg. from pref.)')
% ylabel('response (normalized)')
% title('tuning curve')
% 
% subplot(3,2,2);hold on
% nP = size(score,1);
% clrs = hsv(nP);
% for i = 1:nP
%     pt(i) = plot3(score(i,1),score(i,2),score(i,3),'.','Color',clrs(i,:),'MarkerSize',20);
% end
% plot3([score(:,1);score(1,1)],[score(:,2);score(1,2)],[score(:,3);score(1,3)],'k--')
% xlabel('PC1')
% ylabel('PC2')
% zlabel('PC3')
% legend(pt,num2str(C'))
% title('PCA')
% 
% subplot(3,2,4);hold on
% stem(cumsum(explained))
% xlabel('PC')
% ylabel('cum.sum explained variance')
% ylim([0 100])
% 
% subplot(3,2,5);hold on
% imagesc(D)
% axis tight
% xticks([1:size(D,2)])
% xticklabels(num2str(C'))
% yticks([1:size(D,1)])
% yticklabels(num2str(C'))
% xBox = [0 1 1 0]; xOrthOff = [(0.5:11.5) (0.5:3.5)]; xBox = repmat(xBox',1,length(xOrthOff))+xOrthOff;
% yBox = [0 0 1 1]; yOrthOff = [(4.5:15.5) (12.5:15.5)]; yBox = repmat(yBox',1,length(yOrthOff))+yOrthOff;
% patch(yBox,xBox,'r','FaceColor','none','EdgeColor','r','LineWidth',2)
% 
% subplot(3,2,6)
% plot(C,score(:,1:4),'LineWidth',2)
% % polarplot(deg2rad(C),score(:,1:4),'LineWidth',2)
% % rlim([-1.5 1.5])
% legend({'PC1','PC2','PC3','PC4'})
% 
% sgtitle(ttl)
% saveas(gcf,fullfile(path,figFileName))

