clear all
close all

% dataFold = '/Volumes/NielsenHome2/Brandon/data';
dataFold = 'Y:\Brandon\data';
ageRange = [37 52];

anaMode = 'MU';

switch anaMode
case 'SU'
    load(fullfile(dataFold,'dataSets','cooling','V1cool_ori',anaMode,['V1cool_ori_' anaMode 'dataSet.mat']))
    data{1} = vertcat(dat.cntrl{:,1});
    data{2} = vertcat(dat.cool{:,1});
case 'MU'
    load(fullfile(dataFold,'dataSets','cooling','V1cool_ori',anaMode,['V1cool_ori_' anaMode 'dataSet.mat']))
    data{1} = vertcat(dat.cntrl{:,1});
    data{2} = vertcat(dat.cool{:,1});
case 'matchedSU'
    load(fullfile(dataFold,'dataSets','cooling','V1cool_ori','matchedSU','V1cool_ori_matchedSUdataSet.mat'))
    data{1} = dat{1};
    data{2} = dat{2};
end
lbl{1} = 'control';
eClr{1} = 'k';
lbl{2} = 'cool V1';
eClr{2} = 'c';

dir = load(fullfile(dataFold,'dataSets','DSdev','anaRSA_dir.mat'));
ori = load(fullfile(dataFold,'dataSets','DSdev','anaRSA_ori.mat'));

for e = 1:2
    for u = 1:height(data{e})
        uAge{e}(u) = ages(strcmp(animals,data{e}.exptName{u}(1:5)));
    end
    goodId{e} = screenUnits(data{e},anaMode);
end

for e = 1:2
    switch anaMode
    case 'SU'
        idx = goodId{e} & uAge{e}(:)>=ageRange(1) & uAge{e}(:)<=ageRange(2);
    case 'MU'
        idx = goodId{e} & uAge{e}(:)>=ageRange(1) & uAge{e}(:)<=ageRange(2);
    case 'matchedSU'
        idx = goodId{1}&goodId{2} & uAge{e}(:)>=ageRange(1) & uAge{e}(:)<=ageRange(2);
    end
    nU(e) = sum(idx);
    dist{e} = anaPCA(data{e}(idx,:));
    rMat{e} = compute_rMat(data{e}(idx,:));


    %LDA
    ar = [];
    for b = 1:100
        uIdx = randi(nU(e),1,min(nU,[],'all'));

        r = rMat{e}.rTrial_norm(:,uIdx);
        cDir = rMat{e}.cTrial;
        cOri = mod(cDir,180);
        nObs = length(cDir);

        [ldaAcc_dir_boot(e,b),ldaAccSem_dir_boot(e,b)] = lda_bn(r,cDir);
        [ldaAcc_ori_boot(e,b),ldaAccSem_ori_boot(e,b)] = lda_bn(r,cOri);
        
        [ldaAcc_dir_bootShuff(e,b),ldaAccSem_dir_bootShuff(e,b)] = lda_bn(r,cDir(randperm(nObs,nObs)));
        [ldaAcc_ori_bootShuff(e,b),ldaAccSem_ori_bootShuff(e,b)] = lda_bn(r,cOri(randperm(nObs,nObs)));
    end
    accDir(e) = mean(ldaAcc_dir_boot(e,:));
    accDir_sem(e) = mean(ldaAccSem_dir_boot(e,:));
    accOri(e) = mean(ldaAcc_ori_boot(e,:));
    accOri_sem(e) = mean(ldaAccSem_ori_boot(e,:));

    accDir_shuff(e) = mean(ldaAcc_dir_bootShuff(e,:));
    accDir_sem_shuff(e) = mean(ldaAccSem_dir_bootShuff(e,:));
    accOri_shuff(e) = mean(ldaAcc_ori_bootShuff(e,:));
    accOri_sem_shuff(e) = mean(ldaAccSem_ori_bootShuff(e,:));

%     [accDir(e),accDir_sem(e)] = lda_bn(rMat{e}.rTrial_norm,rMat{e}.cTrial);
%     [accOri(e),accOri_sem(e)] = lda_bn(rMat{e}.rTrial_norm,mod(rMat{e}.cTrial,180));


    %RSA
    %compute the spearman correlation between pairwise distance vectors for
    %template and empirical data
    dir_Corr(e) = corr(dir.diss',dist{e}.diss','type','Spearman');
    dirNull_1{e} = corr(dir.diss',dist{e}.dissNull,'type','Spearman');
    ori_Corr(e) = corr(ori.diss',dist{e}.diss','type','Spearman');
    oriNull_1{e} = corr(ori.diss',dist{e}.dissNull,'type','Spearman');

    %compute the partial spearman correlation between distance vectors
    dir_pCorr(e) = partialcorr(dir.diss',dist{e}.diss',ori.diss','type','Spearman');
    dirNull{e} = partialcorr(dir.diss',dist{e}.dissNull,ori.diss','type','Spearman');
    ori_pCorr(e) = partialcorr(ori.diss',dist{e}.diss',dir.diss','type','Spearman');
    oriNull{e} = partialcorr(ori.diss',dist{e}.dissNull,dir.diss','type','Spearman');

    %estimate the standard deviation of correlation by bootstraping
    %shuffled data
    nBoots = 100;
    nCond = size(dist{e}.rMean,1);
    boots = randi(nCond,[nBoots,nCond]);
    for b = 1:size(boots,1)
        dissBoot = pdist(dist{e}.rMean(boots(b,:),:),'spearman');
        dirBoots_Corr(b) = corr(dir.diss',dissBoot','type','Spearman');
        oriBoots_Corr(b) = corr(ori.diss',dissBoot','type','Spearman');
        dirBoots_pCorr(b) = partialcorr(dir.diss',dissBoot',ori.diss','type','Spearman');
        oriBoots_pCorr(b) = partialcorr(ori.diss',dissBoot',dir.diss','type','Spearman');
    end
    dirStd_pCorr(e) = std(dirBoots_pCorr);
    oriStd_pCorr(e) = std(oriBoots_pCorr);
    dirStd_Corr(e) = std(dirBoots_Corr);
    oriStd_Corr(e) = std(oriBoots_Corr);
    clear dirBoots_pCorr oriBoots_pCorr dirBoots_Corr oriBoots_Corr


end

figure('Position',[100 100 800 500]);
nCol = 3;
for e = 1:2

    subplot(2,nCol,1+(nCol*(e-1)));hold on
    imagesc(rMat{e}.rMean_norm)
    axis square tight
    box on
    title(['n = ' num2str(nU(e))])

    subplot(2,nCol,2+(nCol*(e-1)));hold on
    [coeff,score{e},~,~,explained{e}] = pca(rMat{e}.rMean_norm);
    nC = size(rMat{e}.rMean_norm,1);
    clrs = hsv(nC);
    for i = 1:nC
        plot3(score{e}(i,1),score{e}(i,2),score{e}(i,3),'o','MarkerEdgeColor','none','MarkerFaceColor',clrs(i,:))
    end
    plot3([score{e}(:,1);score{e}(1,1)],[score{e}(:,2);score{e}(1,2)],[score{e}(:,3);score{e}(1,3)],'k')
    view(3)
    axis square
    box on
    title('PCA')

%     subplot(2,nCol,3+(nCol*(e-1)));hold on
%     [coeff,score,~,~,~] = pca(rMat{e}.rTrial_norm);
%     for i = 1:nC
%        idx = rMat{e}.cTrial == rMat{e}.cMean(i);
%        plot3(score(idx,1),score(idx,2),score(idx,3),'o','MarkerEdgeColor','none','MarkerFaceColor',clrs(i,:))
%     end
%     % [score] = tsne(dist{e}.rMean_norm,'NumDimensions',3);
%     % for i = 1:nC
%     %     plot3(score(i,1),score(i,2),score(i,3),'o','MarkerEdgeColor','none','MarkerFaceColor',clrs(i,:))
%     % end
%     % plot3(score(:,1),score(:,2),score(:,3),'k')
%     view(3)
%     axis square
%     box on


    subplot(2,nCol,3+(nCol*(e-1)));hold on
    rdm = squareform(pdist(rMat{e}.rMean_norm,'spearman'));
    imagesc(rdm)
    axis square tight
    box on
    title('RDM')

end
sgtitle(['V1 cool: P' num2str(ageRange(1)) '-' num2str(ageRange(2))])

figure;
nCol = 2;
for e = 1:2


    subplot(2,nCol,1);hold on
    stem(explained{e},'Color',eClr{e})
    xlabel('PC')
    ylabel('variance explained')
    axis square
    box on

    subplot(2,nCol,2);hold on
    errorbar(accDir(e),accOri(e),accOri_sem(e),accOri_sem(e),accDir_sem(e),accDir_sem(e),'Color',eClr{e})
    errorbar(accDir_shuff(e),accOri_shuff(e),accOri_sem_shuff(e),accOri_sem_shuff(e),accDir_sem_shuff(e),accDir_sem_shuff(e),'Color',eClr{e})
    plot([0 1],[0 1],'k--')
    xline(1/12,'k:')
    yline(1/6,'k:')
    title('LDA performace')
    xlabel('dir acc')
    ylabel('ori acc')
    axis square
    box on

    subplot(2,nCol,3);hold on
    plot(1:length(explained{e}),cumsum(explained{e}),'Color',eClr{e})
    xlabel('PC')
    ylabel('cumulative variance explained')
    axis square
    box on

    subplot(2,nCol,4);hold on
    errorbar(dir_Corr(e),ori_Corr(e),oriStd_Corr(e),oriStd_Corr(e),dirStd_Corr(e),dirStd_Corr(e),'Color',eClr{e})
    plot([0 1],[0 1],'k--')
    title('RDM correlations')
    xlabel('dir corr')
    ylabel('ori corr')
    axis square
    box on


end
sgtitle(['V1 cool: P' num2str(ageRange(1)) '-' num2str(ageRange(2))])

