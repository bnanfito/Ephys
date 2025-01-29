



% clear all
% close all
% load('/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/cooling/V1cool_MU_ori/V1cool_SU_ori_projectTbl.mat')
% ageGroup = {[28  30],[31  33],[34 37],[38  40],[46 200]};
% % ageGroup = {[0 33],[34 39],[40 300]};
% % ageGroup = {[28 30],[29 31],[30 32],[31 33],[32 34],[33 35],[34 36],[35 37],[36 38]};
% % ageGroup = {[35 37],[36 38],[37 39],[38 40],[39 41],[40 42],[41 43],[42 44],[43 45]};
% % ageGroup = {[28 30],[31 33],[34 36],[37 39],[40 200]};
% nAG = length(ageGroup);
% animals = unique(projTbl.experimentId);
% nAn = length(animals);
% for a = 1:nAn
%     ages(a) = unique(projTbl.age(strcmp(projTbl.experimentId,animals(a))));
% end
% coolBit = 0;
% coolInd = projTbl.duringMFlag == coolBit;
% 
% figure;
% for a = 1:nAG
%     ageInd = projTbl.age>=ageGroup{a}(1) & projTbl.age<=ageGroup{a}(2);
%     sumStats = vertcat(projTbl.sumStats{ageInd & coolInd});
%     if isempty(sumStats)
%         continue
%     end
%     sumStats = sumStats(sumStats.goodUnit,:);
%     if isempty(sumStats)
%         continue
%     end
%     animals_AG = unique({projTbl.experimentId{ageInd & coolInd}});
%     nA(a) = length(animals_AG);
%     nU(a) = height(sumStats);
% 
%     [x,y,rCent,tCent,score,D,Dshift,distF,distNull] = anaPCA(sumStats,0);
% %     ttl = ['P' num2str(ageGroup{ag}(1)) '-' num2str(ageGroup{ag}(2)) '; Na=' num2str(nA(ag)) '; Nu=' num2str(nU(ag))];
% %     sgtitle(ttl)
% 
% 
% 
%     subplot(5,nAG,a+(nAG*0));hold on
%     imagesc(x);
%     axis tight
%     yticks(1:size(x,1))
%     yticklabels(num2str(y'))
%     colorbar
%     xlabel('unit')
%     ylabel('direction of motion')
% 
%     ttl = ['P' num2str(ageGroup{a}(1)) '-' num2str(ageGroup{a}(2)) '; Na=' num2str(nA(a)) '; Nu=' num2str(nU(a))];
%     title(ttl)
% 
% 
%     subplot(5,nAG,a+(nAG*1));hold on
%     plot(tCent,rCent) 
%     plot(tCent,mean(rCent,2,'omitnan'),'k','LineWidth',2)
%     plot(repmat(tCent,2,1),mean(rCent,2,'omitnan')'+([1;-1].*(std(rCent,[],2,'omitnan')/sqrt(size(rCent,2)))'),'k','LineWidth',2)
%     xlabel('deg. relative to preferred')
%     ylabel('mean (+/-sem) normalized response')
% 
% 
%     subplot(5,nAG,a+(nAG*2));hold on
%     np = size(x,1);
%     clrs = hsv(np);
%     for i = 1:np
%     clr = clrs(i,:);
%     pt(i) = plot3(score(i,1),score(i,2),score(i,3),'.','Color',clr,'MarkerSize',20);
%     end
%     plot3([score(:,1);score(1,1)],[score(:,2);score(1,2)],[score(:,3);score(1,3)],'k--','LineWidth',2)
%     if a == 1
%         legend(pt,num2str(y'))
%     end
%     xlabel('PC1');ylabel('PC2');zlabel('PC3')
% 
% 
%     subplot(5,nAG,a+(nAG*3));hold on
%     imagesc(D);
%     axis tight
%     yticks(1:length(y))
%     yticklabels(num2str(y'))
%     xticks(1:length(y))
%     xticklabels(num2str(y'))
%     colorbar
%     xlabel('direction of motion')
%     ylabel('direction of motion')
% 
% 
%     subplot(5,nAG,a+(nAG*4));hold on
%     angDisp = y(y<=180);
%     plot(angDisp,distF,'k-o','LineWidth',2)
%     sem = std(Dshift)/sqrt(size(Dshift,1)); sem = sem(y<=180);
%     v = var(Dshift); v = v(y<=180);
%     patch([angDisp fliplr(angDisp)],[distF-v fliplr(distF+v)],'k','EdgeColor','none','FaceAlpha',0.2)
%     
%     plot(angDisp,mean(distNull),'r-o')
%     sem = std(distNull)/sqrt(size(distNull,1));
%     v = var(distNull);
%     sig2 = std(distNull,'omitnan')*2;
%     for i = 1:size(distNull,2)
%         P99(i) = prctile(distNull(:,i),99,'Method','exact');
%         P95(i) = prctile(distNull(:,i),95,'Method','exact');
%         P05(i) = prctile(distNull(:,i),5,'Method','exact');
%         P01(i) = prctile(distNull(:,i),1,'Method','exact');
%     end
%     patch([angDisp fliplr(angDisp)],[mean(distNull)-v fliplr(mean(distNull)+v)],'r','EdgeColor','none','FaceAlpha',0.2)
%     patch([angDisp fliplr(angDisp)],[P05 fliplr(P95)],'r','FaceColor','none','EdgeColor','r','LineStyle','--')
%     patch([angDisp fliplr(angDisp)],[P01 fliplr(P99)],'r','FaceColor','none','EdgeColor','r','LineStyle',':')
%         
%     sigHiX = angDisp(distF>P95);
%     sigHiY = distF(distF>P95);
%     sigLoX = angDisp(distF<P05);
%     sigLoY = distF(distF<P05);
%     text(sigHiX,sigHiY+(sigHiY*0.1),'*')
%     text(sigLoX,sigLoY-(sigLoY*0.1),'*')
%     
%     xticks([0 90 180])
%     xticklabels({'0','90','180'})
%     xlabel('angular disparity (+/- deg)')
% 
% 
% 
% end





















%anaV1Cool_ori
clear
close all

anaMode = 'MU';
proj = ['V1cool_MU_ori'];
area = 'PSS';

% dataFold = fullfile('Y:\Brandon\data');
% projTbl = getProjectFiles(proj,1,'age','recSite','penNr','priorMFlag','priorDescr',...
%                                        'duringMFlag','manipDescr','manipDetail',...
%                                        'looperNameCond1','looperNameCond2',...
%                                        'looperNameCond3','looperNameCond4',...
%                                        'looperNameCond5');
% projTbl = projTbl(strcmp(projTbl.recSite,area),:);
% animals = unique(projTbl.experimentId);
% for a = 1:length(animals)
%     aniIdx = strcmp(projTbl.experimentId,animals{a});
%     ages(a) = unique(projTbl.age(aniIdx));
% end
% for e = 1:height(projTbl)
%     animal = projTbl.experimentId{e};
%     unit = projTbl.unitNr{e};
%     expt = projTbl.experimentNr{e};
%     probe = projTbl.probeId(e);
%     exptName = projTbl.fileBase{e};
%     disp(['generating sumStats for ' exptName])
%     [sumStats{e,1}] = anaOri(animal,unit,expt,probe,anaMode,dataFold,0,0);
%     nGoodUnits(e,1) = sum(sumStats{e}.goodUnit);
% end
% projTbl.sumStats = sumStats;
% projTbl.nGU = nGoodUnits;

load(['Y:\Brandon\data\dataSets\cooling\' proj '\V1cool_' anaMode '_ori_projectTbl.mat'])
% load(['/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/cooling/V1cool_MU_ori/V1cool_' anaMode '_ori_projectTbl.mat'])


coolIdx = projTbl.duringMFlag==1 & strcmp(projTbl.manipDescr,'Cooling') & ...
          strcmp(projTbl.manipDetail,'V1');
datTbl = []; %project table without repeated experiments in the same penetration
for a = 1:length(animals)
    aniIdx = strcmp(projTbl.experimentId,animals{a});

    %cntrl
    idx = ~coolIdx & aniIdx;
    cntrlPens = unique(projTbl(idx,:).penNr);
    for pen = 1:length(cntrlPens)
        curDat = projTbl(idx & projTbl.penNr==cntrlPens(pen),:);
        maxGUidx = find( curDat.nGU == max(curDat.nGU) );
        if length(maxGUidx)>1
            maxGUidx = maxGUidx(1);
        end
        datTbl = vertcat(datTbl,curDat(maxGUidx,:));
    end
    

    %cool
    idx = coolIdx & aniIdx;
    coolPens = unique(projTbl(idx,:).penNr);
    for pen = 1:length(coolPens)
        curDat = projTbl(idx & projTbl.penNr==coolPens(pen),:);
        maxGUidx = find( curDat.nGU == max(curDat.nGU) );
        if length(maxGUidx)>1
            maxGUidx = maxGUidx(1);
        end
        datTbl = vertcat(datTbl,curDat(maxGUidx,:));
    end

end

%Organize sumStats column from datTbl by animal and condition in vars
%'cntrlAniDat' and 'coolAniDat'
coolIdx = datTbl.duringMFlag==1 & strcmp(datTbl.manipDescr,'Cooling') & ...
          strcmp(datTbl.manipDetail,'V1');
for a = 1:length(animals)
    aniIdx = strcmp(datTbl.experimentId,animals{a});

    idx = ~coolIdx & aniIdx;
    cntrlAniDat{a} = vertcat(datTbl.sumStats{idx});
    if ~isempty(cntrlAniDat{a})
        cntrlAniDat{a} = cntrlAniDat{a}(cntrlAniDat{a}.goodUnit,:);
    end

    idx = coolIdx & aniIdx;
    coolAniDat{a} = vertcat(datTbl.sumStats{idx});
    if ~isempty(coolAniDat{a})
        coolAniDat{a} = coolAniDat{a}(coolAniDat{a}.goodUnit,:);
    end

end

%Organize sumStats from datTbl by age as defined by var 'ageGroups'
% ageGroups = {[29 32],[33 36],[37 52],[53 300]};
ageGroups = {[29 32],[33 37],[38 300]};
for ag = 1:length(ageGroups)
    agIdx = datTbl.age>=ageGroups{ag}(1) & datTbl.age<=ageGroups{ag}(2);

    idx = ~coolIdx & agIdx;
    cntrlAgeDat{ag} = vertcat(datTbl.sumStats{idx});
    if ~isempty(cntrlAgeDat{ag})
        cntrlAgeDat{ag} = cntrlAgeDat{ag}(cntrlAgeDat{ag}.goodUnit,:);
    end

    idx = coolIdx & agIdx;
    coolAgeDat{ag} = vertcat(datTbl.sumStats{idx});
    if ~isempty(coolAgeDat{ag})
        coolAgeDat{ag} = coolAgeDat{ag}(coolAgeDat{ag}.goodUnit,:);
    end

end


%% Plot

figure;
subplot(2,2,1)
histogram(ages,[20:1:100])
ylabel('animal count')
xlabel('age')
ax1 = gca;

sz1 = 5;
sz2 = 20;
lw = 1;

rPref = nan(length(animals),2);
metric = {'ldr','rPref','latency'};
for m = 1:length(metric)

    subplot(2,2,m+1);hold on
    for a = 1:length(animals)

        x = ages(a);
        for c = 1:2

            if c == 1
                dat = cntrlAniDat{a};  
                clr = 'k';
            elseif c == 2
                dat = coolAniDat{a};
                clr = 'c';
            end
            if ~isempty(dat)
                if strcmp(metric{m},'dsi')
                    Y = dat.dsi;
                elseif strcmp(metric{m},'ldr')
                    Y = dat.ldr;
                elseif strcmp(metric{m},'rPref')
                    Y = dat.rPref;
                    rPref(a,c) = mean(Y,'omitnan');
                elseif strcmp(metric{m},'rNull')
                    Y = dat.rNull;
                elseif strcmp(metric{m},'latency')
                    Y = dat.latency;
                end
                y = mean(Y,'omitnan');
                n = length(Y);
                sem = std(Y,'omitnan')/sqrt(n);
                if n < 5
                    plot(x,y,[clr 'o'],'MarkerSize',sz1,'LineWidth',lw)
                else
                    plot(x,y,[clr '.'],'MarkerSize',sz2)
                end
                plot(repmat(x,2,1),y+([1;-1]*sem),[clr '-'],'LineWidth',lw)
            end

        end
    
    end
    if strcmp(metric{m},'dsi')
        ylabel('DSI')
    elseif strcmp(metric{m},'ldr')
        ylabel('Ldir')
    elseif strcmp(metric{m},'rPref')
        ylabel('response to pref. (Hz)')
    elseif strcmp(metric{m},'rNull')
        ylabel('response to null (Hz)')
    elseif strcmp(metric{m},'latency')
        ylabel('latency')
    end
    xlabel('age')
    ax{m} = gca;

end

linkaxes([ax1 ax{:}],'x')
sgtitle([proj ' ' area ' ' anaMode ' ind. animal means'])

SI = (rPref(:,2) - rPref(:,1)) ./ rPref(:,1);
figure;plot(ages,SI,'o')



%% Distance analysis

figure;
nAG = length(ageGroups);
for a = 1:nAG

    sumStats = coolAgeDat{a};
    if isempty(sumStats)
        continue
    end
    sumStats = sumStats(sumStats.goodUnit,:);
    if isempty(sumStats)
        continue
    end

    names = vertcat(sumStats.exptName{:}); names = unique(names(:,1:5),'rows');
    nA(a) = size(names,1);
    nU(a) = height(sumStats);

    [y,~,x,~,rCent,tCent,score,coeff,D,Dshift,distF,distNull] = anaPCA(sumStats); x = x./max(x);
%     [x,y,rCent,tCent,score,D,Dshift,distF,distNull] = anaPCA(sumStats,0);
%     ttl = ['P' num2str(ageGroup{ag}(1)) '-' num2str(ageGroup{ag}(2)) '; Na=' num2str(nA(ag)) '; Nu=' num2str(nU(ag))];
%     sgtitle(ttl)



    subplot(5,nAG,a+(nAG*0));hold on
    imagesc(x);
    axis tight
    yticks(1:size(x,1))
    yticklabels(num2str(y'))
    colorbar
    xlabel('unit')
    ylabel('direction of motion')

    ttl = ['P' num2str(ageGroups{a}(1)) '-' num2str(ageGroups{a}(2)) '; Na=' num2str(nA(a)) '; Nu=' num2str(nU(a))];
    title(ttl)


    subplot(5,nAG,a+(nAG*1));hold on
    plot(tCent,rCent) 
    plot(tCent,mean(rCent,2,'omitnan'),'k','LineWidth',2)
    plot(repmat(tCent,2,1),mean(rCent,2,'omitnan')'+([1;-1].*(std(rCent,[],2,'omitnan')/sqrt(size(rCent,2)))'),'k','LineWidth',2)
    xlabel('deg. relative to preferred')
    ylabel('mean (+/-sem) normalized response')


    subplot(5,nAG,a+(nAG*2));hold on
    np = size(x,1);
    clrs = hsv(np);
    for i = 1:np
    clr = clrs(i,:);
    pt(i) = plot3(score(i,1),score(i,2),score(i,3),'.','Color',clr,'MarkerSize',20);
    end
    plot3([score(:,1);score(1,1)],[score(:,2);score(1,2)],[score(:,3);score(1,3)],'k--','LineWidth',2)
    if a == 1
        legend(pt,num2str(y'))
    end
    xlabel('PC1');ylabel('PC2');zlabel('PC3')


    subplot(5,nAG,a+(nAG*3));hold on
    imagesc(D);
    axis tight
    yticks(1:length(y))
    yticklabels(num2str(y'))
    xticks(1:length(y))
    xticklabels(num2str(y'))
    colorbar
    xlabel('direction of motion')
    ylabel('direction of motion')


    subplot(5,nAG,a+(nAG*4));hold on
    angDisp = y(y<=180);
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
        
    sigHiX = angDisp(distF>P95);
    sigHiY = distF(distF>P95);
    sigLoX = angDisp(distF<P05);
    sigLoY = distF(distF<P05);
    text(sigHiX,sigHiY+(sigHiY*0.1),'*')
    text(sigLoX,sigLoY-(sigLoY*0.1),'*')
    
    xticks([0 90 180])
    xticklabels({'0','90','180'})
    xlabel('angular disparity (+/- deg)')

%     clear x y rCent tCent score D Dshift distF distNull

end





























% figure;
% nAG = length(ageGroups);
% for ag = 1:nAG
% 
%     dat = cntrlAgeDat{ag};
%     [~,sortIdx] = sort(dat.oriPref);
%     dat = dat(sortIdx,:);
% 
%     eNames = vertcat(dat.exptName{:});
%     nA(ag) = size(unique(eNames(:,1:5),'rows'),1);
%     
%     R{ag} = cat(3,dat.response{:});
%     r{ag} = squeeze(mean(R{ag},1,'omitnan'));
%     r{ag}(r{ag}<0)=0;
%     r{ag} = r{ag}./max(r{ag});
%     c = [0:22.5:337.5];
% 
%     [coeff{ag},score{ag},latent{ag},tsquare{ag}] = pca(r{ag});
%     D{ag} = squareform(pdist(r{ag},'euclidean'));
% 
%     for i = 1:size(D{ag},1)
%         shift = size(D{ag},1)-(i-1);
%         Dshift{ag}(i,:) = circshift(D{ag}(i,:),shift,2);
%     end
% 
%     subplot(4,nAG,ag+(nAG*0))
%     imagesc(r{ag})
%     yticks(1:4:length(c))
%     yticklabels(num2str(c(1:4:end)'))
%     ylabel('direction')
%     xlabel('unit')
%     title(['ages ' num2str(ageGroups{ag}(1)) '-' num2str(ageGroups{ag}(2)) '; nA=' num2str(nA(ag))])
% 
%     subplot(4,nAG,ag+(nAG*1));hold on
%     np = size(score{ag},1);
%     clrs = hsv(np);
%     for i = 1:np
%         pt(i) = plot3(score{ag}(i,1),score{ag}(i,2),score{ag}(i,3),'.','MarkerSize',20,'Color',clrs(i,:));
%     end
%     plot3([score{ag}(:,1);score{ag}(1,1)],[score{ag}(:,2);score{ag}(1,2)],[score{ag}(:,3);score{ag}(1,3)],'k--')
%     xlabel('PC1');ylabel('PC2');zlabel('PC3')
%     if ag == 1
%         legend(pt,num2str(c'))
%     end
% 
%     subplot(4,nAG,ag+(nAG*2))
%     imagesc(D{ag})
%     xticks(1:4:length(c))
%     xticklabels(num2str(c(1:4:end)'))
%     xlabel('dirction')
%     yticks(1:4:length(c))
%     yticklabels(num2str(c(1:4:end)'))
%     ylabel('direction')
% 
%     subplot(4,nAG,ag+(nAG*3))
%     plot(mean(Dshift{ag}))
%     xticks(1:4:length(c))
%     xticklabels(num2str(c(1:4:end)'))
%     xlabel('angle from direction')
% 
% 
% end
% sgtitle([proj ' ' area ' ' anaMode ' cntrl distance'])
% 
% figure;
% nAG = length(ageGroups);
% for ag = 1:nAG
% 
%     dat = coolAgeDat{ag};
%     [~,sortIdx] = sort(dat.oriPref);
%     dat = dat(sortIdx,:);
% 
%     eNames = vertcat(dat.exptName{:});
%     nA(ag) = size(unique(eNames(:,1:5),'rows'),1);
%     
%     R{ag} = cat(3,dat.response{:});
%     r{ag} = squeeze(mean(R{ag},1,'omitnan'));
%     r{ag}(r{ag}<0)=0;
%     r{ag} = r{ag}./max(r{ag});
%     c = [0:22.5:337.5];
% 
%     [coeff{ag},score{ag},latent{ag},tsquare{ag}] = pca(r{ag});
%     D{ag} = squareform(pdist(r{ag},'euclidean'));
% 
%     for i = 1:size(D{ag},1)
%         shift = size(D{ag},1)-(i-1);
%         Dshift{ag}(i,:) = circshift(D{ag}(i,:),shift,2);
%     end
% 
%     subplot(4,nAG,ag+(nAG*0))
%     imagesc(r{ag})
%     yticks(1:4:length(c))
%     yticklabels(num2str(c(1:4:end)'))
%     ylabel('direction')
%     xlabel('unit')
%     title(['ages ' num2str(ageGroups{ag}(1)) '-' num2str(ageGroups{ag}(2)) '; nA=' num2str(nA(ag))])
% 
%     subplot(4,nAG,ag+(nAG*1));hold on
%     np = size(score{ag},1);
%     clrs = hsv(np);
%     for i = 1:np
%         pt(i) = plot3(score{ag}(i,1),score{ag}(i,2),score{ag}(i,3),'.','MarkerSize',20,'Color',clrs(i,:));
%     end
%     plot3([score{ag}(:,1);score{ag}(1,1)],[score{ag}(:,2);score{ag}(1,2)],[score{ag}(:,3);score{ag}(1,3)],'k--')
%     xlabel('PC1');ylabel('PC2');zlabel('PC3')
% 
%     subplot(4,nAG,ag+(nAG*2))
%     imagesc(D{ag})
%     xticks(1:4:length(c))
%     xticklabels(num2str(c(1:4:end)'))
%     xlabel('dirction')
%     yticks(1:4:length(c))
%     yticklabels(num2str(c(1:4:end)'))
%     ylabel('direction')
% 
%     subplot(4,nAG,ag+(nAG*3))
%     plot(mean(Dshift{ag}))
%     xticks(1:4:length(c))
%     xticklabels(num2str(c(1:4:end)'))
%     ylabel('mean distance (eucl.)')
%     xlabel('angle bw d1 & d2')
% 
% 
% end
% sgtitle([proj ' ' area ' ' anaMode ' cooled distance'])
% 
% %pref ori dist.
% figure;
% for ag = 1:nAG
%     subplot(nAG,1,ag);hold on
% 
%     dat = cntrlAgeDat{ag};
%     dat = dat(dat.ldr>0.1,:);
%     for u = 1:height(dat)
% 
%         oriPref(u) = dat.meanVec{u}.angDir;
% 
%     end
%     histogram(oriPref,c)
%     clear oriPref
% 
% end




