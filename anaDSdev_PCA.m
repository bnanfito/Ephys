clear all
close all

anaMode = 'SU';
dataFold = '/Volumes/Lab drive/Brandon/data/dataSets/DSdev';
load(fullfile(dataFold,"DSdev_projectTbl.mat"))

ageGroups = {[0 30],[31 33],[34 36],[37 39],[40 300]};
nAG = length(ageGroups);
for ag = 1:nAG
    ageLims = ageGroups{ag};
    v1Dat{ag} = vertcat(projectTbl(strcmp(projectTbl.recSite,'V1') & projectTbl.age>=ageLims(1) & projectTbl.age<=ageLims(2),:).sumStats{:});
    pssDat{ag} = vertcat(projectTbl(strcmp(projectTbl.recSite,'PSS') & projectTbl.age>=ageLims(1) & projectTbl.age<=ageLims(2),:).sumStats{:});

    v1Dat{ag} = v1Dat{ag}(v1Dat{ag}.goodUnit,:);
    pssDat{ag} = pssDat{ag}(pssDat{ag}.goodUnit,:);

    nU_v1(ag) = height(v1Dat{ag});
    for u = 1:nU_v1(ag)
        r = mean(v1Dat{ag}.response{u},'omitnan');
        c = v1Dat{ag}.condition{u}(strcmp(v1Dat{ag}.paramKey{u},'ori'),:);
        if length(r)~=12
            cInt = 0:30:330;
            rInt = interp1([c c(1)+360],[r r(1)],cInt);
            c = cInt;
            r = rInt;
        end
        rV1{ag}(:,u) = r./max(r);
        cV1{ag}(:,u) = c;
    end

    nU_pss(ag) = height(pssDat{ag});
    for u = 1:nU_pss(ag)
        r = mean(pssDat{ag}.response{u},'omitnan');
        c = pssDat{ag}.condition{u}(strcmp(pssDat{ag}.paramKey{u},'ori'),:);
        if length(r)~=12
            cInt = 0:30:330;
            rInt = interp1([c c(1)+360],[r r(1)],cInt);
            c = cInt;
            r = rInt;
        end
        rPSS{ag}(:,u) = r./max(r);
        cPSS{ag}(:,u) = c;
    end

    [coeff{1,ag}, score{1,ag}, latent{1,ag}, tsq{1,ag}, explained{1,ag}] = pca(rV1{ag});
    [coeff{2,ag}, score{2,ag}, latent{2,ag}, tsq{2,ag}, explained{2,ag}] = pca(rPSS{ag});

end


figure;hold on
for ag = 1:nAG
    subplot(2,nAG,ag);hold on
    nP = size(score{1,ag},1);
    clrs = hsv(nP);
    for i = 1:nP
        pt(i) = plot3(score{1,ag}(i,1),score{1,ag}(i,2),score{1,ag}(i,3),'.','Color',clrs(i,:),'MarkerSize',20);
    end
    plot3(score{1,ag}(:,1),score{1,ag}(:,2),score{1,ag}(:,3),'k--')
    title(['age:' num2str(ageGroups{ag}(1)) '-' num2str(ageGroups{ag}(2)) ';nU=' num2str(nU_v1(ag))])
    if ag ==nAG
        legend(pt,num2str(c'))
    end

    subplot(2,nAG,ag+nAG);hold on
    nP = size(score{2,ag},1);
    clrs = hsv(nP);
    for i = 1:nP
        plot3(score{2,ag}(i,1),score{2,ag}(i,2),score{2,ag}(i,3),'.','Color',clrs(i,:),'MarkerSize',20)
    end
    plot3(score{2,ag}(:,1),score{2,ag}(:,2),score{2,ag}(:,3),'k--')
    title(['age:' num2str(ageGroups{ag}(1)) '-' num2str(ageGroups{ag}(2)) ';nU=' num2str(nU_pss(ag))])

end




