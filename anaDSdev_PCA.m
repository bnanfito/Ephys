clear all
% close all

anaMode = 'SU';
% dataFold = '/Volumes/Lab drive/Brandon/data/dataSets/DSdev';
% dataFold = 'F:\Brandon\data\dataSets\DSdev';
dataFold = 'Y:\Brandon\data\dataSets\DSdev';
load(fullfile(dataFold,"DSdev_projectTbl.mat"))

area = 'V1';
ageGroups = {[0 33],[34 39],[40 300]};
nAG = length(ageGroups);
for ag = 1:nAG
    ageLims = ageGroups{ag};
    areaIdx = strcmp(projectTbl.recSite,area);
    ageLimIdx = projectTbl.age>=ageLims(1) & projectTbl.age<=ageLims(2);
    dat{ag} = vertcat(projectTbl(areaIdx & ageLimIdx,:).sumStats{:});

    dat{ag} = dat{ag}(dat{ag}.goodUnit,:);

    [~,sortIdx] = sort(dat{ag}.oriPref);
    dat{ag} = dat{ag}(sortIdx,:);

    nU(ag) = height(dat{ag});
    for u = 1:nU(ag)
        rTmp = mean(dat{ag}.response{u},'omitnan');
        c = dat{ag}.condition{u}(strcmp(dat{ag}.paramKey{u},'ori'),:);
        if length(rTmp)~=12
            cInt = 0:30:330;
            rInt = interp1([c c(1)+360],[rTmp rTmp(1)],cInt);
            c = cInt;
            rTmp = rInt;
        end
        rTmp = rTmp./max(rTmp);
        rTmp(rTmp<0) = 0;
        r{ag}(:,u) = rTmp;
    end

    [coeff{1,ag}, score{1,ag}, latent{1,ag}, tsq{1,ag}, explained{1,ag}] = pca(r{ag});

%     D{1,ag} = squareform(pdist(r{ag}));
    D{1,ag} = dist(r{ag}');

    for i = 1:length(c)
        
        angDif{ag}(i,:) = c-c(i);
        angDif{ag}(i,angDif{ag}(i,:)>180) = angDif{ag}(i,angDif{ag}(i,:)>180)-360;
        angDif{ag}(i,angDif{ag}(i,:)<-180) = angDif{ag}(i,angDif{ag}(i,:)<-180)+360;
        
        a = r{ag}(i,:);
        for j = 1:length(angDif{ag})
            b = r{ag}(j,:);
            D2{1,ag}(i,j) = norm(a-b);
        end
    end

    for i = 1:size(D{ag},1)
        shift = size(D{ag},1)-(i-1);
        Dshift{ag}(i,:) = circshift(D{ag}(i,:),shift,2);
    end

    for i = 1:size(D2{ag},1)
        shift = size(D2{ag},1)-(i-1);
        D2shift{ag}(i,:) = circshift(D2{ag}(i,:),shift,2);
    end

    for i = 1:size(angDif{ag},1)
        shift = size(angDif{ag},1)-(i-1);
        angDifShift{ag}(i,:) = circshift(angDif{ag}(i,:),shift,2);
    end


end


figure;hold on
for ag = 1:nAG

    subplot(4,nAG,ag+(nAG*0));hold on
    title([area ' age:' num2str(ageGroups{ag}(1)) '-' num2str(ageGroups{ag}(2)) ';nU=' num2str(nU(ag))])
    imagesc(r{ag});
    axis tight
    yticks([1:size(r{ag},1)])
    yticklabels(num2str(c'))

    subplot(4,nAG,ag+(nAG*1));hold on
    imagesc(D{1,ag})
    axis tight
    xticks([1:size(D{1,ag},2)])
    xticklabels(num2str(c'))
    yticks([1:size(D{1,ag},1)])
    yticklabels(num2str(c'))

    subplot(4,nAG,ag+(nAG*2));hold on
    nP = size(score{1,ag},1);
    clrs = hsv(nP);
    for i = 1:nP
        pt(i) = plot3(score{1,ag}(i,1),score{1,ag}(i,2),score{1,ag}(i,3),'.','Color',clrs(i,:),'MarkerSize',20);
    end
    plot3([score{1,ag}(:,1);score{1,ag}(1,1)],[score{1,ag}(:,2);score{1,ag}(1,2)],[score{1,ag}(:,3);score{1,ag}(1,3)],'k--')
    xlabel('PC1');ylabel('PC2');zlabel('PC3')
    if ag ==nAG
        legend(pt,num2str(c'))
    end

    subplot(4,nAG,ag+(nAG*3));hold on
    plot(c,mean(Dshift{ag}))
    axis tight
    xticks([0 90 180])
    xlim([0 180])

end




