clear all
close all

%% Load data

anaMode = 'MU';

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
    %sort units in data table by their preferred direction
    [~,sortIdx] = sort(dat{ar,ag}.oriPref);
    dat{ar,ag} = dat{ar,ag}(sortIdx,:);
    nU(ar,ag) = height(dat{ar,ag});
    D = dat{ar,ag};

    rMat{ar,ag} = compute_rMat(D);

    [distDat{ar,ag}] = anaPCA(D);

    %tuning metrics
    ldr(ar,ag) = mean(D.ldr);
    ldr_sem(ar,ag) = sem(D.ldr);
    dsi(ar,ag) = mean(D.dsi);
    dsi_sem(ar,ag) = sem(D.dsi);
    lor(ar,ag) = mean(D.lor);
    lor_sem(ar,ag) = sem(D.lor);
    osi(ar,ag) = mean(D.osi);
    osi_sem(ar,ag) = sem(D.osi);

end
end

%% RSA

for ar = 1:nAR
for ag = 1:nAG

    %compute the spearman correlation between pairwise distance vectors for
    %template and empirical data
    dir_Corr(ar,ag) = corr(dir.diss',distDat{ar,ag}.diss','type','Spearman');
    dirNull_1{ar,ag} = corr(dir.diss',distDat{ar,ag}.dissNull,'type','Spearman');
    ori_Corr(ar,ag) = corr(ori.diss',distDat{ar,ag}.diss','type','Spearman');
    oriNull_1{ar,ag} = corr(ori.diss',distDat{ar,ag}.dissNull,'type','Spearman');

    %compute the partial spearman correlation between distance vectors
    dir_pCorr(ar,ag) = partialcorr(dir.diss',distDat{ar,ag}.diss',ori.diss','type','Spearman');
    dirNull{ar,ag} = partialcorr(dir.diss',distDat{ar,ag}.dissNull,ori.diss','type','Spearman');
    ori_pCorr(ar,ag) = partialcorr(ori.diss',distDat{ar,ag}.diss',dir.diss','type','Spearman');
    oriNull{ar,ag} = partialcorr(ori.diss',distDat{ar,ag}.dissNull,dir.diss','type','Spearman');

    %estimate the standard deviation of correlation by bootstraping
    %shuffled data
    nBoots = 100;
    nCond = size(distDat{ar,ag}.rMean,1);
    boots = randi(nCond,[nBoots,nCond]);
    for b = 1:size(boots,1)
        dissBoot = pdist(distDat{ar,ag}.rMean(boots(b,:),:),'spearman');
        dirBoots_Corr(b) = corr(dir.diss',dissBoot','type','Spearman');
        oriBoots_Corr(b) = corr(ori.diss',dissBoot','type','Spearman');
        dirBoots_pCorr(b) = partialcorr(dir.diss',dissBoot',ori.diss','type','Spearman');
        oriBoots_pCorr(b) = partialcorr(ori.diss',dissBoot',dir.diss','type','Spearman');
    end
    dirStd_pCorr(ar,ag) = std(dirBoots_pCorr);
    oriStd_pCorr(ar,ag) = std(oriBoots_pCorr);
    dirStd_Corr(ar,ag) = std(dirBoots_Corr);
    oriStd_Corr(ar,ag) = std(oriBoots_Corr);
    clear dirBoots_pCorr oriBoots_pCorr dirBoots_Corr oriBoots_Corr

end
end

%% LDA

for ar = 1:length(areas)
    for ag = 1:length(ageGroups)
        for b = 1:100
%                 uIdx = randperm(nU(ar,ag),nSmpl(nS));
            uIdx = randi(nU(ar,ag),1,min(nU,[],'all'));

            r = rMat{ar,ag}.rTrial_norm(:,uIdx);
            cDir = rMat{ar,ag}.cTrial;
            cOri = mod(cDir,180);
            nObs = length(cDir);

            [ldaAcc_dir_boot(ar,ag,b),ldaAccSem_dir_boot(ar,ag,b)] = lda_bn(r,cDir);
            [ldaAcc_ori_boot(ar,ag,b),ldaAccSem_ori_boot(ar,ag,b)] = lda_bn(r,cOri);
            
            [ldaAcc_dir_bootShuff(ar,ag,b),ldaAccSem_dir_bootShuff(ar,ag,b)] = lda_bn(r,cDir(randperm(nObs,nObs)));
            [ldaAcc_ori_bootShuff(ar,ag,b),ldaAccSem_ori_bootShuff(ar,ag,b)] = lda_bn(r,cOri(randperm(nObs,nObs)));
        end
        ldaAcc_dir(ar,ag) = mean(ldaAcc_dir_boot(ar,ag,:));
        ldaAccSem_dir(ar,ag) = mean(ldaAccSem_dir_boot(ar,ag,:));
        ldaAcc_ori(ar,ag) = mean(ldaAcc_ori_boot(ar,ag,:));
        ldaAccSem_ori(ar,ag) = mean(ldaAccSem_ori_boot(ar,ag,:));

        ldaAcc_dir_shuff(ar,ag) = mean(ldaAcc_dir_bootShuff(ar,ag,:));
        ldaAccSem_dir_shuff(ar,ag) = mean(ldaAccSem_dir_bootShuff(ar,ag,:));
        ldaAcc_ori_shuff(ar,ag) = mean(ldaAcc_ori_bootShuff(ar,ag,:));
        ldaAccSem_ori_shuff(ar,ag) = mean(ldaAccSem_ori_bootShuff(ar,ag,:));
    end
end

% figure; hold on
% count = 0;
% linStyl = {'-','--'};
% curAG = 3;
% for ar = 1:nAR
%     if strcmp(areas{ar},'V1')
%         clr = 'b';
%         lblArea = 'V1';
%     elseif strcmp(areas{ar},'PSS')
%         clr = 'r';
%         lblArea = 'PSS';
%     end
% for ag = curAG
% 
% 
%     D = distDat{ar,ag}; R = D.rTrial_norm; 
% %     R = R-mean(R); R = R./max(R);
% %     R = zscore(R);
%     C_dir = D.cTrial; dirs = unique(C_dir);
%     C_ori = mod(C_dir,180); oris = unique(C_ori);
%     [coeff,score] = pca(R);
% 
%     for pc = 1:size(score,2)
%         [acc{ar,ag,1}(pc)] = lda_bn(score(:,1:pc),C_dir);
%         [acc{ar,ag,2}(pc)] = lda_bn(score(:,1:pc),C_ori);
%     end
%     count = count+1;
%     p(count) = plot(acc{ar,ag,1},[clr linStyl{1}],'LineWidth',2);
%     lbl{count} = [lblArea ' dir'];
%     count = count+1;
%     p(count) = plot(acc{ar,ag,2},[clr linStyl{2}],'LineWidth',2);
%     lbl{count} = [lblArea ' ori'];
% end
% end
% title(['P' num2str(ageGroups{curAG}(1)) ' - ' num2str(ageGroups{curAG}(2))])
% legend(p,lbl,'Location','southeast')
% ylabel('accuracy')
% xlabel('number of PCs')
% clear lbl
% 
% 
% D = distDat{2,3}; R = D.rTrial_norm; 
% % R = R./max(R); 
% C_dir = D.cTrial; dirs = unique(C_dir);
% C_ori = mod(C_dir,180); oris = unique(C_ori);
% [coeff,score] = pca(R);
% % score = D.score;
% pcA = 1;pcB = 2;pcC = 3;
% mdl_dir2 = fitcdiscr(score(:,[pcA pcB]),C_dir);
% mdl_dir3 = fitcdiscr(score(:,[pcA pcB pcC]),C_dir);
% mdl_ori2 = fitcdiscr(score(:,[pcA pcB]),C_ori);
% mdl_ori3 = fitcdiscr(score(:,[pcA pcB pcC]),C_ori);
% 
% eval = 90;
% a = find(dirs==eval);
% b = find(dirs~=eval)';
% sym(a) = 'o'; sym(b) = 'o';
% clrs = hsv(length(dirs));
% 
% figure;hold on
% for i = b
% K = mdl_dir2.Coeffs(a,i).Const; L = mdl_dir2.Coeffs(a,i).Linear; f = @(x1,x2) K + L(1)*x1+L(2)*x2;
% h2 = fimplicit(f); h2.LineWidth = 2; h2.Color = clrs(i,:); h2.DisplayName = ['boundary between ' num2str(dirs(a)) '&' num2str(dirs(i))];
% end
% for i = 1:length(dirs)
%     idx = C_dir==dirs(i);
%     p1(i) = plot(score(idx,pcA),score(idx,pcB),sym(i),'MarkerEdgeColor','k','MarkerFaceColor',clrs(i,:));
%     lbl{i} = num2str(dirs(i));
% end
% % gs = gscatter(score(:,pcA),score(:,pcB),C_dir,[],sym);
% xlabel(['PC' num2str(pcA)]);ylabel(['PC' num2str(pcB)])
% legend(p1,lbl)
% xlim([-2 2])
% ylim([-2 2])
% axis square
% box on
% 
% figure;hold on
% for i = b(1)
% K = mdl_dir3.Coeffs(a,i).Const; L = mdl_dir3.Coeffs(a,i).Linear; f = @(x1,x2,x3) K + L(1)*x1+L(2)*x2+L(3)*x3;
% h2 = fimplicit3(f); h2.EdgeColor = 'none'; h2.FaceColor = clrs(i,:); h2.FaceAlpha = 0.3; h2.DisplayName = ['boundary between ' num2str(dirs(a)) '&' num2str(dirs(i))];
% end
% for i = 1:length(dirs)
%     idx = C_dir==dirs(i);
%     p2(i) = plot3(score(idx,pcA),score(idx,pcB),score(idx,pcC),sym(i),'MarkerEdgeColor','k','MarkerFaceColor',clrs(i,:));
% end
% xlabel(['PC' num2str(pcA)]);ylabel(['PC' num2str(pcB)]);zlabel(['PC' num2str(pcC)])
% legend(p2,lbl)
% view([62,8.2])
% xlim([-2 2])
% ylim([-1.5 1.5])
% zlim([-1.5 1.5])
% axis square
% box on
% 
% clear sym lbl
% 
% eval = mod(eval,180);
% a = find(oris==eval);
% b = find(oris~=eval)';
% sym(a) = 'o'; sym(b) = 'o';
% clrs = hsv(length(oris));
% 
% figure;hold on
% for i = b
%     K = mdl_ori2.Coeffs(a,i).Const; L = mdl_ori2.Coeffs(a,i).Linear; f = @(x1,x2) K + L(1)*x1+L(2)*x2;
%     h2 = fimplicit(f); h2.LineWidth = 2; h2.Color = clrs(i,:); h2.DisplayName = ['boundary between ' num2str(oris(a)) '&' num2str(oris(i))];
% end
% for i = 1:length(oris)
%     idx = C_ori==oris(i);
%     p3(i) = plot(score(idx,pcA),score(idx,pcB),sym(i),'MarkerEdgeColor','k','MarkerFaceColor',clrs(i,:));
%     lbl{i} = num2str(oris(i));
% end
% % gs = gscatter(score(:,pcA),score(:,pcB),C_ori,[],sym);
% xlabel(['PC' num2str(pcA)]);ylabel(['PC' num2str(pcB)])
% legend(p3,lbl)
% xlim([-2 2])
% ylim([-2 2])
% axis square
% box on
% 
% figure;hold on
% for i = b(1)
% K = mdl_ori3.Coeffs(a,i).Const; L = mdl_ori3.Coeffs(a,i).Linear; f = @(x1,x2,x3) K + L(1)*x1+L(2)*x2+L(3)*x3;
% h2 = fimplicit3(f); h2.EdgeColor = 'none'; h2.FaceColor = clrs(i,:); h2.FaceAlpha = 0.3; h2.DisplayName = ['boundary between ' num2str(oris(a)) '&' num2str(oris(i))];
% end
% for i = 1:length(oris)
%     idx = C_ori==oris(i);
%     p4(i) = plot3(score(idx,pcA),score(idx,pcB),score(idx,pcC),sym(i),'MarkerEdgeColor','k','MarkerFaceColor',clrs(i,:));
% end
% xlabel(['PC' num2str(pcA)]);ylabel(['PC' num2str(pcB)]);zlabel(['PC' num2str(pcC)])
% legend(p4,lbl)
% view([62,8.2])
% xlim([-2 2])
% ylim([-1.5 1.5])
% zlim([-1.5 1.5])
% axis square
% box on

%% Population vector

%trialNum = 38;
for trialNum = 1:60
ar = 2;ag = 3;
t = rMat{ar,ag}.cPref(2,:);
c = rMat{ar,ag}.cTrial(trialNum);
r = rMat{ar,ag}.rTrial_norm(trialNum,:);
mv = meanvec(t,r);
err(trialNum) = rad2deg(angdiff(deg2rad(mv.angDir),deg2rad(c)));
figure
polarplot( deg2rad(t) , r ,'o');hold on
polarplot(repmat(deg2rad(mv.angDir),2,1),[0;max(r)],'k','LineWidth',2)
polarplot(repmat(deg2rad(c),2,1),[0;max(r)],'r','LineWidth',2)
title(num2str(err(trialNum)))
end
figure;polarhistogram(deg2rad(err))

%% Plot

close all
for ar = 1:length(areas)
for ag = 1:length(ageGroups)
    
    figure(1);
    subplot(nAR,nAG,ag+(nAG*(ar-1)));
    imagesc(rMat{ar,ag}.rTrial_norm);
    colormap gray
    axis square
    box on
    
    figure(2);
    subplot(nAR,nAG,ag+(nAG*(ar-1)));
    for b = 1:100
        uIdx = randi(height(dat{ar,ag}),1,min(nU,[],'all')); %bootstrap sample with replacement
        rdm(:,:,b) = squareform(pdist(rMat{ar,ag}.rMean_norm(:,uIdx),'spearman'));
    end
    imagesc(mean(rdm,3));
    colormap gray
    axis square
    box on
    
    figure(3);
    subplot(nAR,nAG,ag+(nAG*(ar-1)));hold on
    [coeff,score] = pca(rMat{ar,ag}.rTrial_norm);
    clrs = hsv(length(rMat{ar,ag}.cMean));
    pc_A = 1;pc_B = 2;pc_C = 4;
    if size(score,1)==size(rMat{ar,ag}.rTrial,1)
        for c = 1:size(clrs,1)
        idx = rMat{ar,ag}.cTrial==rMat{ar,ag}.cMean(c);
        plot3(score(idx,pc_A),score(idx,pc_B),score(idx,pc_C),'o','Color',clrs(c,:),'MarkerFaceColor',clrs(c,:))
        end
    else
        plot3([score(:,pc_A);score(1,pc_A)],[score(:,pc_B);score(1,pc_B)],[score(:,pc_C);score(1,pc_C)],'k')
        for c = 1:size(clrs,1)
        plot3(score(c,pc_A),score(c,pc_B),score(c,pc_C),'o','Color',clrs(c,:),'MarkerFaceColor',clrs(c,:))
        end
    end
    view(3)
    axis square
    box on
end
end


% RSA

f2 = figure;
f2.Position = [100 100 800 800];

subplot(2,2,1)
imagesc(ori.rdm)
colormap gray

subplot(2,2,2); hold on
b1 = bar(1:numel(ori_pCorr),ori_pCorr(:));
b1.FaceColor = 'flat';
v1Idx = find(strcmp(areas,'V1')):2:numel(ori_pCorr);
pssIdx = find(strcmp(areas,'PSS')):2:numel(ori_pCorr);
b1.CData(v1Idx,:) = repmat([0 0 1],length(v1Idx),1);
b1.CData(pssIdx,:) = repmat([1 0 0],length(pssIdx),1);
ylabel('partial Spearman corr.')
xticks([1.5:2:numel(ori_pCorr)])
for ag = 1:nAG
    xtLbl{ag} = ['P' num2str(ageGroups{ag}(1)) '-' num2str(ageGroups{ag}(2))];
end
xticklabels(xtLbl)
er1 = errorbar(1:numel(ori_pCorr),ori_pCorr(:),oriStd_pCorr(:),oriStd_pCorr(:));
er1.Color = [0 0 0];
er1.LineStyle = 'none';

subplot(2,2,3)
imagesc(dir.rdm)
colormap gray

subplot(2,2,4); hold on
b2 = bar(1:numel(dir_pCorr),dir_pCorr(:));
b2.FaceColor = 'flat';
b2.CData(v1Idx,:) = repmat([0 0 1],length(v1Idx),1);
b2.CData(pssIdx,:) = repmat([1 0 0],length(pssIdx),1);
ylabel('partial Spearman corr.')
xticks([1.5:2:numel(dir_pCorr)])
for ag = 1:nAG
    xtLbl{ag} = ['P' num2str(ageGroups{ag}(1)) '-' num2str(ageGroups{ag}(2))];
end
xticklabels(xtLbl)
er2 = errorbar(1:numel(dir_pCorr),dir_pCorr(:),dirStd_pCorr(:),dirStd_pCorr(:));
er2.Color = [0 0 0];
er2.LineStyle = 'none';


figure;hold on
for a = 1:length(areas)
if strcmp(areas{a},'V1')
clr = 'b';
elseif strcmp(areas{a},'PSS')
clr = 'r';
end
errorbar(ori_Corr(a,:),dir_Corr(a,:),dirStd_Corr(a,:),dirStd_Corr(a,:),oriStd_Corr(a,:),oriStd_Corr(a,:),[clr '--'])
errorbar(ori_pCorr(a,:),dir_pCorr(a,:),dirStd_pCorr(a,:),dirStd_pCorr(a,:),oriStd_pCorr(a,:),oriStd_pCorr(a,:),[clr 'o-'],'LineWidth',2)
end
xlabel('spearman corr. with orientation')
ylabel('spearman corr. with direction')
plot([-0.4 1],[-0.4 1],'k--')
axis square
box on
legend({'V1','V1 (partial)','PMLS','PMLS (partial)'},'Location','southwest')

figure;hold on
subplot(2,2,1);hold on
plot([0 1],[0 1],'k--')
subplot(2,2,2);hold on
plot([0 1],[0 1],'k--')
subplot(2,2,3);hold on
plot([-0.5 1],[-0.5 1],'k--')
subplot(2,2,4);hold on
plot([-0.2 1],[-0.2 1],'k--')
for ar = 1:length(areas)
    switch areas{ar}
        case 'V1'
            clr = 'b';
        case 'PSS'
            clr = 'r';
    end
    for ag = 1:length(ageGroups)
        D = dat{ar,ag};
    end
    
    subplot(2,2,1);hold on
    p(ar) = errorbar(ldr(ar,:),lor(ar,:),lor_sem(ar,:),lor_sem(ar,:),ldr_sem(ar,:),ldr_sem(ar,:),[clr]);
    xlabel('LDR'); xlim([0 0.6])
    ylabel('LOR'); ylim([0 0.6])
    axis square
    box on

    subplot(2,2,2);hold on
    errorbar(dsi(ar,:),osi(ar,:),osi_sem(ar,:),osi_sem(ar,:),dsi_sem(ar,:),dsi_sem(ar,:),[clr]);
    xlabel('DSI'); xlim([0.2 1])
    ylabel('OSI'); ylim([0.2 1])
    axis square
    box on

    subplot(2,2,3);hold on
    errorbar(dir_Corr(ar,:),ori_Corr(ar,:),oriStd_Corr(ar,:),oriStd_Corr(ar,:),dirStd_Corr(ar,:),dirStd_Corr(ar,:),[clr])
    xlabel('dir. corr.'); xlim([-0.5 1])
    ylabel('ori. corr.'); ylim([-0.5 1])
    axis square
    box on

    subplot(2,2,4);hold on
    errorbar(dir_pCorr(ar,:),ori_pCorr(ar,:),oriStd_pCorr(ar,:),oriStd_pCorr(ar,:),dirStd_pCorr(ar,:),dirStd_pCorr(ar,:),[clr])
    xlabel('dir. pcorr.'); xlim([-0.2 1])
    ylabel('ori. pcorr.'); ylim([-0.2 1])
    axis square
    box on
end
subplot(2,2,1);hold on
legend(p(1:2),areas,'Location','southeast')


figure; hold on
plot([0 1],[0 1],'k--')
xline(1/12,'k:')
yline(1/6,'k:')
for ar = 1:length(areas)
    switch areas{ar}
        case 'V1'
            clr = 'b';
        case 'PSS'
            clr = 'r';
    end
    errorbar(ldaAcc_dir(ar,:),ldaAcc_ori(ar,:),ldaAccSem_ori(ar,:),ldaAccSem_ori(ar,:),ldaAccSem_dir(ar,:),ldaAccSem_dir(ar,:),clr)
    errorbar(ldaAcc_dir_shuff(ar,:),ldaAcc_ori_shuff(ar,:),ldaAccSem_ori_shuff(ar,:),ldaAccSem_ori_shuff(ar,:),ldaAccSem_dir_shuff(ar,:),ldaAccSem_dir_shuff(ar,:),[clr '--'])
end
axis square
box on
xlabel('LDA acc. on dir.')
xlim([0 1])
ylabel('LDA acc. on ori.')
ylim([0 1])



