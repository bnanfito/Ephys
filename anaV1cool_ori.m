%anaV1Cool_ori
clear
close all

anaMode = 'MU';
proj = ['V1cool_' anaMode '_ori'];
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

% load('Y:\Brandon\data\dataSets\cooling\V1cool_MU_ori\V1cool_MU_ori_projectTbl.mat')
load('/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/cooling/V1cool_MU_ori/V1cool_MU_ori_projectTbl.mat')


coolIdx = projTbl.duringMFlag==1 & strcmp(projTbl.manipDescr,'Cooling') & ...
          strcmp(projTbl.manipDetail,'V1');
datTbl = [];
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

ageGroups = {[28 29],[30 31],[32 33],[34 35],[36 37],[38 45],[46 200]};
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

subplot(2,2,2);hold on
for a = 1:length(animals)
    x = ages(a);
    
    if ~isempty(cntrlAniDat{a})
        Y = cntrlAniDat{a}.ldr;
        y = mean(Y);
        n = length(Y);
        sem = std(Y)/sqrt(n);
        plot(x,y,'ko','MarkerSize',n/5,'LineWidth',2)
        plot(repmat(x,2,1),y+([1;-1]*sem),'k-','LineWidth',2)
    end

    if ~isempty(coolAniDat{a})
        Y = coolAniDat{a}.ldr;
        y = mean(Y);
        n = length(Y);
        sem = std(Y)/sqrt(n);
        plot(x,y,'co','MarkerSize',n/5,'LineWidth',2)
        plot(repmat(x,2,1),y+([1;-1]*sem),'c-','LineWidth',2)
    end

end
ylabel('Ldir')
xlabel('age')
ax2 = gca;

subplot(2,2,3);hold on
for a = 1:length(animals)
    x = ages(a);
    
    if ~isempty(cntrlAniDat{a})
        Y = cntrlAniDat{a}.rNull;
        y = mean(Y);
        n = length(Y);
        sem = std(Y)/sqrt(n);
        plot(x,y,'ko','MarkerSize',n/5,'LineWidth',2)
        plot(repmat(x,2,1),y+([1;-1]*sem),'k-','LineWidth',2)
    end

    if ~isempty(coolAniDat{a})
        Y = coolAniDat{a}.rNull;
        y = mean(Y);
        n = length(Y);
        sem = std(Y)/sqrt(n);
        plot(x,y,'co','MarkerSize',n/5,'LineWidth',2)
        plot(repmat(x,2,1),y+([1;-1]*sem),'c-','LineWidth',2)
    end

end
ylabel('response to null Hz')
xlabel('age')
ax3 = gca;

subplot(2,2,4);hold on
for a = 1:length(animals)
    x = ages(a);
    
    if ~isempty(cntrlAniDat{a})
        Y = cntrlAniDat{a}.rPref;
        y = mean(Y);
        n = length(Y);
        sem = std(Y)/sqrt(n);
        plot(x,y,'ko','MarkerSize',n/5,'LineWidth',2)
        plot(repmat(x,2,1),y+([1;-1]*sem),'k-','LineWidth',2)
    end

    if ~isempty(coolAniDat{a})
        Y = coolAniDat{a}.rPref;
        y = mean(Y);
        n = length(Y);
        sem = std(Y)/sqrt(n);
        plot(x,y,'co','MarkerSize',n/5,'LineWidth',2)
        plot(repmat(x,2,1),y+([1;-1]*sem),'c-','LineWidth',2)
    end

end
ylabel('response to pref (Hz)')
xlabel('age')
ax4 = gca;

linkaxes([ax1 ax2 ax3 ax4],'x')
sgtitle([proj ' ' area ' ' anaMode ' ind. animal means'])




figure;
nAG = length(ageGroups);
for ag = 1:nAG

    dat = cntrlAgeDat{ag};
    [~,sortIdx] = sort(dat.oriPref);
    dat = dat(sortIdx,:);

    eNames = vertcat(dat.exptName{:});
    nA(ag) = size(unique(eNames(:,1:5),'rows'),1);
    
    R{ag} = cat(3,dat.response{:});
    r{ag} = squeeze(mean(R{ag},1,'omitnan'));
    r{ag}(r{ag}<0)=0;
    r{ag} = r{ag}./max(r{ag});
    c = [0:22.5:337.5];

    [coeff{ag},score{ag},latent{ag},tsquare{ag}] = pca(r{ag});
    D{ag} = squareform(pdist(r{ag},'euclidean'));

    for i = 1:size(D{ag},1)
        shift = size(D{ag},1)-(i-1);
        Dshift{ag}(i,:) = circshift(D{ag}(i,:),shift,2);
    end

    subplot(4,nAG,ag+(nAG*0))
    imagesc(r{ag})
    yticks(1:4:length(c))
    yticklabels(num2str(c(1:4:end)'))
    ylabel('direction')
    xlabel('unit')
    title(['ages ' num2str(ageGroups{ag}(1)) '-' num2str(ageGroups{ag}(2)) '; nA=' num2str(nA(ag))])

    subplot(4,nAG,ag+(nAG*1))
    plot3(score{ag}(:,1),score{ag}(:,2),score{ag}(:,3),'ko-')
    xlabel('PC1');ylabel('PC2');zlabel('PC3')

    subplot(4,nAG,ag+(nAG*2))
    imagesc(D{ag})
    xticks(1:4:length(c))
    xticklabels(num2str(c(1:4:end)'))
    xlabel('dirction')
    yticks(1:4:length(c))
    yticklabels(num2str(c(1:4:end)'))
    ylabel('direction')

    subplot(4,nAG,ag+(nAG*3))
    plot(mean(Dshift{ag}))
    xticks(1:4:length(c))
    xticklabels(num2str(c(1:4:end)'))
    xlabel('angle from direction')


end
sgtitle([proj ' ' area ' ' anaMode ' cntrl distance'])



figure;
nAG = length(ageGroups);
for ag = 1:nAG

    dat = coolAgeDat{ag};
    [~,sortIdx] = sort(dat.oriPref);
    dat = dat(sortIdx,:);

    eNames = vertcat(dat.exptName{:});
    nA(ag) = size(unique(eNames(:,1:5),'rows'),1);
    
    R{ag} = cat(3,dat.response{:});
    r{ag} = squeeze(mean(R{ag},1,'omitnan'));
    r{ag}(r{ag}<0)=0;
    r{ag} = r{ag}./max(r{ag});
    c = [0:22.5:337.5];

    [coeff{ag},score{ag},latent{ag},tsquare{ag}] = pca(r{ag});
    D{ag} = squareform(pdist(r{ag},'euclidean'));

    for i = 1:size(D{ag},1)
        shift = size(D{ag},1)-(i-1);
        Dshift{ag}(i,:) = circshift(D{ag}(i,:),shift,2);
    end

    subplot(4,nAG,ag+(nAG*0))
    imagesc(r{ag})
    yticks(1:4:length(c))
    yticklabels(num2str(c(1:4:end)'))
    ylabel('direction')
    xlabel('unit')
    title(['ages ' num2str(ageGroups{ag}(1)) '-' num2str(ageGroups{ag}(2)) '; nA=' num2str(nA(ag))])

    subplot(4,nAG,ag+(nAG*1))
    plot3(score{ag}(:,1),score{ag}(:,2),score{ag}(:,3),'ko-')
    xlabel('PC1');ylabel('PC2');zlabel('PC3')

    subplot(4,nAG,ag+(nAG*2))
    imagesc(D{ag})
    xticks(1:4:length(c))
    xticklabels(num2str(c(1:4:end)'))
    xlabel('dirction')
    yticks(1:4:length(c))
    yticklabels(num2str(c(1:4:end)'))
    ylabel('direction')

    subplot(4,nAG,ag+(nAG*3))
    plot(mean(Dshift{ag}))
    xticks(1:4:length(c))
    xticklabels(num2str(c(1:4:end)'))
    ylabel('mean distance (eucl.)')
    xlabel('angle bw d1 & d2')


end
sgtitle([proj ' ' area ' ' anaMode ' cooled distance'])



