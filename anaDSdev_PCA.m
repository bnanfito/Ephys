clear all
close all

anaMode = 'SU';
% dataFold = '/Volumes/Lab drive/Brandon/data/dataSets/DSdev';
% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/DSdev';
dataFold = 'F:\Brandon\data\dataSets\DSdev';
load(fullfile(dataFold,"DSdev_projectTbl.mat"))
% dataFold = 'Y:\Brandon\data\dataSets\DSdev';
% load(fullfile(dataFold,"DSdev_dataSet.mat"))

area = 'PSS';
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

    [r{ag},c,~,~,score{ag},D{1,ag},Dshift{1,ag},distF{1,ag},distNull{1,ag}] = anaPCA(dat{ag},0);

%     for u = 1:nU(ag)
%         rTmp = mean(dat{ag}.response{u},'omitnan');
%         c = dat{ag}.condition{u}(strcmp(dat{ag}.paramKey{u},'ori'),:);
%         if length(rTmp)~=12
%             cInt = 0:30:330;
%             rInt = interp1([c c(1)+360],[rTmp rTmp(1)],cInt);
%             c = cInt;
%             rTmp = rInt;
%         end
%         rTmp = rTmp./max(rTmp);
%         rTmp(rTmp<0) = 0;
%         r{ag}(:,u) = rTmp;
%     end
% 
%     [coeff{1,ag}, score{1,ag}, latent{1,ag}, tsq{1,ag}, explained{1,ag}] = pca(r{ag});
% 
% %     D{1,ag} = squareform(pdist(r{ag}));
%     D{1,ag} = dist(r{ag}');
% 
%     for i = 1:length(c)
%         
%         angDif{ag}(i,:) = c-c(i);
%         angDif{ag}(i,angDif{ag}(i,:)>180) = angDif{ag}(i,angDif{ag}(i,:)>180)-360;
%         angDif{ag}(i,angDif{ag}(i,:)<-180) = angDif{ag}(i,angDif{ag}(i,:)<-180)+360;
%         
%         a = r{ag}(i,:);
%         for j = 1:length(angDif{ag})
%             b = r{ag}(j,:);
%             D2{1,ag}(i,j) = norm(a-b);
%         end
%     end
% 
%     for i = 1:size(D{ag},1)
%         shift = size(D{ag},1)-(i-1);
%         Dshift{ag}(i,:) = circshift(D{ag}(i,:),shift,2);
%     end
% 
%     for i = 1:size(D2{ag},1)
%         shift = size(D2{ag},1)-(i-1);
%         D2shift{ag}(i,:) = circshift(D2{ag}(i,:),shift,2);
%     end
% 
%     for i = 1:size(angDif{ag},1)
%         shift = size(angDif{ag},1)-(i-1);
%         angDifShift{ag}(i,:) = circshift(angDif{ag}(i,:),shift,2);
%     end


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
    imagesc(D{ag})
    axis tight
    xticks([1:size(D{ag},2)])
    xticklabels(num2str(c'))
    yticks([1:size(D{ag},1)])
    yticklabels(num2str(c'))

    subplot(4,nAG,ag+(nAG*2));hold on
    nP = size(score{ag},1);
    clrs = hsv(nP);
    for i = 1:nP
        pt(i) = plot3(score{ag}(i,1),score{ag}(i,2),score{ag}(i,3),'.','Color',clrs(i,:),'MarkerSize',20);
    end
    plot3([score{ag}(:,1);score{ag}(1,1)],[score{ag}(:,2);score{ag}(1,2)],[score{ag}(:,3);score{ag}(1,3)],'k--')
    xlabel('PC1');ylabel('PC2');zlabel('PC3')
    if ag ==nAG
        legend(pt,num2str(c'))
    end

    subplot(4,nAG,ag+(nAG*3));hold on

    angDisp = c(c<=180);
    plot(angDisp,distF{ag},'k-o','LineWidth',2)
    sem = std(Dshift{ag})/sqrt(size(Dshift{ag},1)); sem = sem(c<=180);
    v = var(Dshift{ag}); v = v(c<=180);
    patch([angDisp fliplr(angDisp)],[distF{ag}-v fliplr(distF{ag}+v)],'k','EdgeColor','none','FaceAlpha',0.2)

    plot(angDisp,mean(distNull{ag}),'r-o')
    sem = std(distNull{ag})/sqrt(size(distNull{ag},1));
    v = var(distNull{ag});
    sig2 = std(distNull{ag},'omitnan')*2;
    for i = 1:size(distNull{ag},2)
        P99(i) = prctile(distNull{ag}(:,i),99,'Method','exact');
        P95(i) = prctile(distNull{ag}(:,i),95,'Method','exact');
        P05(i) = prctile(distNull{ag}(:,i),5,'Method','exact');
        P01(i) = prctile(distNull{ag}(:,i),1,'Method','exact');
    end
    patch([angDisp fliplr(angDisp)],[mean(distNull{ag})-v fliplr(mean(distNull{ag})+v)],'r','EdgeColor','none','FaceAlpha',0.2)
    patch([angDisp fliplr(angDisp)],[P05 fliplr(P95)],'r','FaceColor','none','EdgeColor','r','LineStyle','--')
    patch([angDisp fliplr(angDisp)],[P01 fliplr(P99)],'r','FaceColor','none','EdgeColor','r','LineStyle',':')

    sigHiX = angDisp(distF{ag}>P95);
    sigHiY = distF{ag}(distF{ag}>P95);
    sigLoX = angDisp(distF{ag}<P05);
    sigLoY = distF{ag}(distF{ag}<P05);
    text(sigHiX,sigHiY+(sigHiY*0.1),'*')
    text(sigLoX,sigLoY-(sigLoY*0.1),'*')
    
    xticks([0 90 180])
    xticklabels({'0','90','180'})
    xlabel('angular disparity (+/- deg)')

%     plot(c,mean(Dshift{ag}))
%     axis tight
%     xticks([0 90 180])
%     xlim([0 180])

end




