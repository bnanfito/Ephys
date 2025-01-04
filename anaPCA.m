% close all
% clear

function [x,y,rCent,tCent,score,D,Dshift,distF,distNull] = anaPCA(sumStats,plt)

tAve = 0;
fullDim = 1;

% animal = 'febn2';
% unit = '000';
% expt = '002';
% area = 'V1';
% anaMode = 'SU';
% % dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
% % dataFold = '/Volumes/Lab drive/Brandon/data';
% dataFold = '/Volumes/NielsenHome2/Brandon/data';
% % dataFold = 'Y:\Brandon\data';
% 
% [sumStats] = anaOri(animal,unit,expt,area,anaMode,dataFold,0,0);

% load('/Volumes/NielsenHome2/Brandon/data/dataSets/cooling/V1cool_MU_ori/V1cool_MU_ori_projectTbl.mat')
% sumStats = vertcat(projTbl.sumStats{projTbl.age>=30 & projTbl.age<=32 & projTbl.duringMFlag == 0});

sumStats = sumStats(sumStats.goodUnit,:);
[~,oriPrefIdx] = sort(sumStats.oriPref);
sumStats = sumStats(oriPrefIdx,:);

nU = height(sumStats);
for u = 1:nU
    C = sumStats.condition{u}(strcmp(sumStats.paramKey{1},'ori'),:);
    if length(C)==12
        rTrial(:,:,u) = sumStats.response{u}(1:5,:);
    else
        cInt = 0:30:330;
        rTmp = sumStats.response{u}(1:5,:);
        rInt = interp1([C C(1)+360],[rTmp rTmp(:,1)]',cInt);
        rTrial(:,:,u) = rInt';
        C = cInt;
    end
    rTrial(rTrial<0)=0;
    rMean(:,u) = mean(rTrial(:,:,u),'omitnan');
end

nReps = size(rTrial,1);
nConds = size(rTrial,2);
rTrial = reshape(rTrial,nReps*nConds,nU);
% rMean = squeeze(mean(R,1,'omitnan'));
c = repmat(C,nReps,1);c = c(:)';

if tAve == 1
    x = rMean;
    y = C;
else
    x = rTrial;
    y = c;
end
np = size(x,1);
if tAve == 1
    clrs = hsv(np);
else
    clrsTmp = hsv(nConds);
    clrs = [];
    for i = 1:nConds
        clrs = vertcat(clrs,repmat(clrsTmp(i,:),nReps,1));
    end
end
x = x./max(x);

[coeff,score,latent,tsquare,explained] = pca(x);
if fullDim == 1
    D = dist(x');
    % D = squareform(pdist(x,'squaredeuclidean'));
else
    D = dist(score(:,1:3)');
end

for u = 1:nU
    [tCent,rCent(:,u)] = alignDirTuning(y,x(:,u)');
end

%% distance function

for i = 1:size(D,1)
    shift = size(D,1)-(i-1);
    Dshift(i,:) = circshift(D(i,:),shift,2);
end
distF = mean(Dshift);
distF = distF(y<=180);

%% Null Dist

nNullRep = 1000;
for nr = 1:nNullRep
    
%     randIdx = randperm(size(rTrial,1));
%     rTshuff = rTrial(randIdx,:);
%     rMshuff = squeeze(mean(reshape(rTshuff,nReps,nConds,nU),1,'omitnan'));
%     if tAve == 1
%         shuff = rMshuff;
%     else
%         shuff = rTshuff;
%     end

    if tAve == 1
        randIdx = randperm(size(rMean,1));
        rMshuff = rMean(randIdx,:);
        shuff = rMshuff;
    else
        randIdx = randperm(size(rTrial,1));
        rTshuff = rTrial(randIdx,:);
        shuff = rTshuff;
    end

    shuff = shuff./max(shuff);
    [coeffShuff,scoreShuff,latentShuff,tsquareShuff,explainedShuff] = pca(shuff);
    if fullDim == 1
        Dnull(:,:,nr) = dist(shuff');
    else
        Dnull(:,:,nr) = dist(scoreShuff(:,1:3)');
    end
    for i = 1:size(Dnull,1)
        shift = size(Dnull,1)-(i-1);
        DnullShift(i,:,nr) = circshift(Dnull(i,:,nr),shift,2);
    end
    distNull2{nr} = DnullShift(:,:,nr);
    distNull(nr,:) = mean(DnullShift(:,:,nr),1,'omitnan');

end
distNull2 = vertcat(distNull2{:}); distNull2 = distNull2(:,y<=180);
distNull = distNull(:,y<=180);

distF_z = (distF-mean(distNull,'omitnan'))./std(distNull,'omitnan');

%% Plot

if plt == 1

figure('Position',[0 0 1000 1500]);

subplot(4,2,1);hold on;
imagesc(x);
axis tight
if tAve == 1
    yticks(1:size(x,1))
    yticklabels(num2str(y'))
else
    yticks(1:nReps:size(x,1))
    yticklabels(num2str(y(1:nReps:size(x,1))'))
end
colorbar
xlabel('unit')
ylabel('direction of motion')

subplot(4,2,8);hold on;
plot(tCent,rCent) 
plot(tCent,mean(rCent,2,'omitnan'),'k','LineWidth',2)
plot(repmat(tCent,2,1),mean(rCent,2,'omitnan')'+([1;-1].*(std(rCent,[],2,'omitnan')/sqrt(size(rCent,2)))'),'k','LineWidth',2)
xlabel('deg. relative to preferred')
ylabel('mean (+/-sem) normalized response')


subplot(4,2,3);hold on;
imagesc(D);
axis tight
yticks(1:length(y))
yticklabels(num2str(y'))
xticks(1:length(y))
xticklabels(num2str(y'))
colorbar
xBox = [0 1 1 0]; xOrthOff = [(0.5:11.5) (0.5:3.5)]; xBox = repmat(xBox',1,length(xOrthOff))+xOrthOff;
yBox = [0 0 1 1]; yOrthOff = [(4.5:15.5) (12.5:15.5)]; yBox = repmat(yBox',1,length(yOrthOff))+yOrthOff;
patch(yBox,xBox,'r','FaceColor','none','EdgeColor','r','LineWidth',2)
xlabel('direction of motion')
ylabel('direction of motion')
title('RDM')


subplot(4,2,2);hold on;
% for i = 1:np
% clr = clrs(i,:);
% pt(i) = plot3(scoreShuff(i,1),scoreShuff(i,2),scoreShuff(i,3),'o','Color',clr,'MarkerSize',10,'LineWidth',1.5);
% end
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

subplot(4,2,4);hold on
stem(explained)
xlabel('PC')
ylabel('explained variance')
% ylim([0 100])

subplot(4,2,6)
if tAve == 1
    plot(y,score(:,1:4),'LineWidth',2)
else
    plot(y,score(:,1:4),'o')
end
% polarplot(deg2rad(C),score(:,1:4),'LineWidth',2)
% rlim([-1.5 1.5])
legend({'PC1','PC2','PC3','PC4'})

subplot(4,2,7);hold on
angDisp = y(y<=180);
% plot(angDisp,distF_z,'k-o','LineWidth',2)
plot(angDisp,distF,'k-o','LineWidth',2)
sem = std(Dshift)/sqrt(size(Dshift,1)); sem = sem(y<=180);
v = var(Dshift); v = v(y<=180);
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

subplot(4,2,5)
% imagesc(mean(Dnull,3,'omitnan'))
imagesc(Dnull(:,:,end))
axis tight
yticks(1:length(y))
yticklabels(num2str(y'))
xticks(1:length(y))
xticklabels(num2str(y'))
colorbar
xBox = [0 1 1 0]; xOrthOff = [(0.5:11.5) (0.5:3.5)]; xBox = repmat(xBox',1,length(xOrthOff))+xOrthOff;
yBox = [0 0 1 1]; yOrthOff = [(4.5:15.5) (12.5:15.5)]; yBox = repmat(yBox',1,length(yOrthOff))+yOrthOff;
patch(yBox,xBox,'r','FaceColor','none','EdgeColor','r','LineWidth',2)
xlabel('direction of motion')
ylabel('direction of motion')
title('ex. shuffle RDM')

% ttl = [animal ' ' unit ' ' expt ' ' area ' ' anaMode];
% if fullDim == 1
%     ttl = [ttl ' full dim.'];
% else
%     ttl = [ttl ' first 3 PC'];
% end
% sgtitle(ttl)
% figName = [animal '_u' unit '_' expt '_' area '_' anaMode '_distance.fig'];
% saveas(gcf,fullfile(dataFold,'Ephys',animal,[animal '_u' unit '_' expt],figName))

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

