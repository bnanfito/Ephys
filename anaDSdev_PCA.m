clear all
close all

%% Load data

anaMode = 'SU';

% dataFold = '/Volumes/Lab drive/Brandon/data/dataSets/DSdev';
% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/DSdev';
% dataFold = '/Volumes/NielsenHome2/Brandon/data/dataSets/DSdev';
% dataFold = 'F:\Brandon\data\dataSets\DSdev';
dataFold = 'Y:\Brandon\data\dataSets\DSdev';
load(fullfile(dataFold,['DSdev_' anaMode 'dataSet.mat']))
dir = load(fullfile(dataFold,'anaRSA_dir.mat'));
ori = load(fullfile(dataFold,'anaRSA_ori.mat'));

%% Organize data by age/area

areas = {'V1','PSS'};
ageGroups = {[29 32],[33 36],[37 max(projTbl.age)]};
% ageGroups = {[28 32],[29 33],[30 34],[31 35],[32 36],[33 37],[34 38],[35 39],[36 40],[37 41],[38 42],[39 43],[40 44],[41 300]};

nAG = length(ageGroups);
nAR = length(areas);
for ar = 1:nAR
for ag = 1:nAG

    ageLims = ageGroups{ag};
    areaIdx = strcmp(projTbl.recSite,areas{ar});
    ageLimIdx = projTbl.age>=ageLims(1) & projTbl.age<=ageLims(2);
    dat{ar,ag} = vertcat(projTbl(areaIdx & ageLimIdx,:).sumStats{:});

    %only take units that pass inclusion criteria (screenUnits)
    dat{ar,ag} = dat{ar,ag}(screenUnits(dat{ar,ag},anaMode),:);

    [~,sortIdx] = sort(dat{ar,ag}.oriPref);
    dat{ar,ag} = dat{ar,ag}(sortIdx,:);
    nU(ar,ag) = height(dat{ar,ag});

    [distDat{ar,ag}] = anaPCA(dat{ar,ag});

%     [~,guess{ar,ag},truth{ar,ag},confusion{ar,ag}] = popDecode(distDat{ar,ag}.rTrial,distDat{ar,ag}.cTrial);

end
end

%% LDA

% for ar = 1:nAR
% for ag = 1:nAG
% 
%     D = distDat{ar,ag};
%     [acc{ar,ag}] = lda_bn(D.rTrial,D.cTrial);
% 
% end
% end
% 
% figure;hold on
% plot(v1(:,1),'b-')
% plot(v1(:,2),'b--')
% plot(pss(:,1),'r-')
% plot(pss(:,2),'r--')
% yline(1/12,'k-')
% yline(1/6,'k--')
% legend({'v1 dir','v1 ori','pss dir','pss ori','chance dir','chance ori'})
% ylabel('5-fold cv mean acc.')
% xlabel('age group')

D = distDat{1,3}; R = D.rTrial; R = R-mean(R); R = R./max(R); C = D.cTrial;
[coeff,score] = pca(R);
figure;hold on
gscatter(score(:,1),score(:,2),C)
mdl = fitcdiscr(score(:,1:2),C);
a = 9; b = 1:length(unique(C)); b = b(b~=a);
for i = b
K = mdl.Coeffs(a,i).Const; L = mdl.Coeffs(a,i).Linear; f = @(x1,x2) K + L(1)*x1+L(2)*x2;
h2 = fimplicit(f); h2.DisplayName = ['boundary between ' num2str(D.cMean(a)) '&' num2str(D.cMean(i))];
end

%% RSA

for ar = 1:nAR
for ag = 1:nAG

    dirCorr_1(ar,ag) = corr(dir.diss',distDat{ar,ag}.diss','type','Spearman');
    dirNull_1{ar,ag} = corr(dir.diss',distDat{ar,ag}.dissNull,'type','Spearman');
    oriCorr_1(ar,ag) = corr(ori.diss',distDat{ar,ag}.diss','type','Spearman');
    oriNull_1{ar,ag} = corr(ori.diss',distDat{ar,ag}.dissNull,'type','Spearman');

    dirCorr(ar,ag) = partialcorr(dir.diss',distDat{ar,ag}.diss',ori.diss','type','Spearman');
    dirNull{ar,ag} = partialcorr(dir.diss',distDat{ar,ag}.dissNull,ori.diss','type','Spearman');
    oriCorr(ar,ag) = partialcorr(ori.diss',distDat{ar,ag}.diss',dir.diss','type','Spearman');
    oriNull{ar,ag} = partialcorr(ori.diss',distDat{ar,ag}.dissNull,dir.diss','type','Spearman');

    nBoots = 100;
    nCond = size(distDat{ar,ag}.rMean,1);
    boots = randi(nCond,[nBoots,nCond]);
    for b = 1:size(boots,1)
        dissBoot = pdist(distDat{ar,ag}.rMean(boots(b,:),:),'spearman');
        dirBoots_1(b) = corr(dir.diss',dissBoot','type','Spearman');
        oriBoots_1(b) = corr(ori.diss',dissBoot','type','Spearman');
        dirBoots(b) = partialcorr(dir.diss',dissBoot',ori.diss','type','Spearman');
        oriBoots(b) = partialcorr(ori.diss',dissBoot',dir.diss','type','Spearman');
    end
    dirStd(ar,ag) = std(dirBoots);
    oriStd(ar,ag) = std(oriBoots);
    dirStd_1(ar,ag) = std(dirBoots_1);
    oriStd_1(ar,ag) = std(oriBoots_1);
    clear dirBoots oriBoots dirBoots_1 oriBoots_1

end
end


%% Plot

f = figure;
f.Position = [100 100 1000 600];
for ar = 1:nAR
    if strcmp(areas{ar},'V1')
        clr = 'b';
    elseif strcmp(areas{ar},'PSS')
        clr = 'r';
    end
    for ag = 1:nAG
    
        curDat = distDat{ar,ag};
        subplot(nAR,nAG,ag+(nAG*(ar-1)))

        imagesc(curDat.rTrial)
        xlabel('neuron #')
        ylabel('trial #')

        title([areas{ar} '; P' num2str(ageGroups{ag}(1)) '-' num2str(ageGroups{ag}(2))])
    end
end


f = figure;
f.Position = [100 100 1000 600];
for ar = 1:nAR
    if strcmp(areas{ar},'V1')
        clr = 'b';
    elseif strcmp(areas{ar},'PSS')
        clr = 'r';
    end
    for ag = 1:nAG
    
        curDat = distDat{ar,ag};
        subplot(nAR,nAG,ag+(nAG*(ar-1)))

        imagesc(curDat.rdm)
        colorbar
        if ag == 1
            nYtick = length(curDat.cMean);
            yticks(1:2:nYtick)
            yticklabels(curDat.cMean(1:2:nYtick))
        else
            yticks([])
        end
        if ar == nAR
            nXtick = length(curDat.cMean);
            xticks(1:2:nXtick)
            xticklabels(curDat.cMean(1:2:nXtick))
        else
            xticks([])
        end

        title([areas{ar} '; P' num2str(ageGroups{ag}(1)) '-' num2str(ageGroups{ag}(2))])
    end
end


f = figure;
f.Position = [100 100 1000 600];
for ar = 1:nAR
    if strcmp(areas{ar},'V1')
        clr = 'b';
    elseif strcmp(areas{ar},'PSS')
        clr = 'r';
    end
    for ag = 1:nAG
    
        curDat = distDat{ar,ag};
        subplot(nAR,nAG,ag+(nAG*(ar-1)))
        
        
        hold on
        conds = unique(curDat.cMean);
        score = vertcat(curDat.score,curDat.score(1,:));
        cID = curDat.cMean;
%         score = mdscale(curDat.rdm,3); score = vertcat(score,score(1,:));
        plot3(score(:,1),score(:,2),score(:,3),'k','LineWidth',1.5)
        np = size(score);
        dirClrs = hsv(length(conds));
        oriClrs = repmat(hsv(length(conds)/2),2,1);
        for i = conds
            idx = cID == i;
            plot3(score(idx,1),score(idx,2),score(idx,3),'.','Color',dirClrs(conds==i,:),'MarkerSize',30)
        end
        xlabel('PC1')
        ylabel('PC2')
        zlabel('PC3')
        view(3)
        box on
%         grid on
        oriID = mod(curDat.cMean,180)';
        oris = unique(oriID);
%         for o = 1:length(oris)
%             plot3(curDat.score(oriID==oris(o),1),curDat.score(oriID==oris(o),2),curDat.score(oriID==oris(o),3),'--','Color',[0.6 0.6 0.6],'LineWidth',1.5)
%         end

        title([areas{ar} '; P' num2str(ageGroups{ag}(1)) '-' num2str(ageGroups{ag}(2))])
    end
end


f = figure;
f.Position = [100 100 1000 600];
for ar = 1:nAR
    if strcmp(areas{ar},'V1')
        clr = 'b';
    elseif strcmp(areas{ar},'PSS')
        clr = 'r';
    end
    for ag = 1:nAG
    
        curDat = distDat{ar,ag};
        subplot(nAR,nAG,ag+(nAG*(ar-1)))

        hold on
        dispIdx = curDat.cMean>0 & curDat.cMean<=180;
        angDisp = curDat.cMean(dispIdx);
        distF = curDat.distF(dispIdx);
        distNull = curDat.distNull(:,dispIdx);
        plot(angDisp,distF,'k-o','LineWidth',2)
        sem = std(curDat.rdmShift)/sqrt(size(curDat.rdmShift,1)); sem = sem(dispIdx);
        v = var(curDat.rdmShift); v = v(dispIdx);
        patch([angDisp fliplr(angDisp)],[distF-sem fliplr(distF+sem)],'k','EdgeColor','none','FaceAlpha',0.2)
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
        patch([angDisp fliplr(angDisp)],[mean(distNull)-sem fliplr(mean(distNull)+sem)],'r','EdgeColor','none','FaceAlpha',0.2)
        patch([angDisp fliplr(angDisp)],[P05 fliplr(P95)],'r','FaceColor','none','EdgeColor','r','LineStyle','--')
        patch([angDisp fliplr(angDisp)],[P01 fliplr(P99)],'r','FaceColor','none','EdgeColor','r','LineStyle',':')
        sigHiX = angDisp(distF>P95);
        sigHiY = distF(distF>P95);
        sigLoX = angDisp(distF<P05);
        sigLoY = distF(distF<P05);
        text(sigHiX,sigHiY+(sigHiY*0.1),'*')
        text(sigLoX,sigLoY-(sigLoY*0.1),'*')
        xticks([min(angDisp) 90 max(angDisp)])
        xlim([min(angDisp) max(angDisp)])
        if ar == nAR
            xlabel('angular disparity (+/- deg)')
        end
        if ag == 1
            ylabel('Euclidean distance')
        end

%         hold on
%         hDirNull = histogram(dirNull{ar,ag});
%         hDirNull.FaceColor = 'g';
%         xline(dirCorr(ar,ag),'g')
%         xline(prctile(dirNull{ar,ag},95),'g--')
%         xline(prctile(dirNull{ar,ag},99),'g-.')
%         xline(prctile(dirNull{ar,ag},99.9),'g:')
%         hOriNull = histogram(oriNull{ar,ag});
%         hOriNull.FaceColor = 'r';
%         xline(oriCorr(ar,ag),'r')
%         xline(prctile(oriNull{ar,ag},95),'r--')
%         xline(prctile(oriNull{ar,ag},99),'r-.')
%         xline(prctile(oriNull{ar,ag},99.9),'r:')
%         legend([hDirNull hOriNull],{'dir','ori'})

        title([areas{ar} '; P' num2str(ageGroups{ag}(1)) '-' num2str(ageGroups{ag}(2))])
    end
end


% RSA

f2 = figure;
f2.Position = [100 100 800 800];

subplot(2,2,1)
imagesc(ori.rdm)

subplot(2,2,2); hold on
b1 = bar(1:numel(oriCorr),oriCorr(:));
b1.FaceColor = 'flat';
v1Idx = find(strcmp(areas,'V1')):2:numel(oriCorr);
pssIdx = find(strcmp(areas,'PSS')):2:numel(oriCorr);
b1.CData(v1Idx,:) = repmat([0 0 1],length(v1Idx),1);
b1.CData(pssIdx,:) = repmat([1 0 0],length(pssIdx),1);
ylabel('partial Spearman corr.')
xticks([1.5:2:numel(oriCorr)])
for ag = 1:nAG
    xtLbl{ag} = ['P' num2str(ageGroups{ag}(1)) '-' num2str(ageGroups{ag}(2))];
end
xticklabels(xtLbl)
er1 = errorbar(1:numel(oriCorr),oriCorr(:),oriStd(:),oriStd(:));
er1.Color = [0 0 0];
er1.LineStyle = 'none';

subplot(2,2,3)
imagesc(dir.rdm)

subplot(2,2,4); hold on
b2 = bar(1:numel(dirCorr),dirCorr(:));
b2.FaceColor = 'flat';
b2.CData(v1Idx,:) = repmat([0 0 1],length(v1Idx),1);
b2.CData(pssIdx,:) = repmat([1 0 0],length(pssIdx),1);
ylabel('partial Spearman corr.')
xticks([1.5:2:numel(dirCorr)])
for ag = 1:nAG
    xtLbl{ag} = ['P' num2str(ageGroups{ag}(1)) '-' num2str(ageGroups{ag}(2))];
end
xticklabels(xtLbl)
er2 = errorbar(1:numel(dirCorr),dirCorr(:),dirStd(:),dirStd(:));
er2.Color = [0 0 0];
er2.LineStyle = 'none';


figure;hold on
for a = 1:length(areas)
if strcmp(areas{a},'V1')
clr = 'b';
elseif strcmp(areas{a},'PSS')
clr = 'r';
end
errorbar(oriCorr_1(a,:),dirCorr_1(a,:),dirStd_1(a,:),dirStd_1(a,:),oriStd_1(a,:),oriStd_1(a,:),[clr '--'])
errorbar(oriCorr(a,:),dirCorr(a,:),dirStd(a,:),dirStd(a,:),oriStd(a,:),oriStd(a,:),[clr 'o-'],'LineWidth',2)
end
xlabel('spearman corr. with orientation')
ylabel('spearman corr. with direction')
plot([-0.4 1],[-0.4 1],'k--')
legend({'V1','V1 (partial)','PMLS','PMLS (partial)'})
