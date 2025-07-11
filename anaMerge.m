%anaMerge
clear all
close all

proj = 'V1cool_ori';
% dataFold = '/Volumes/NielsenHome2/Brandon/data';
dataFold = 'Y:\Brandon\data';
% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
anaMode = 'SU';

%% Build Project Table

% projTbl = getProjectFiles(proj,1,'age','recSite','penNr','priorMFlag','priorDescr',...
%                                        'duringMFlag','manipDescr','manipDetail',...
%                                        'looperNameCond1','looperNameCond2',...
%                                        'looperNameCond3','looperNameCond4',...
%                                        'looperNameCond5');
% animals = unique(projTbl.experimentId);
% sumStats = [];
% for a = 1:length(animals)
%     aniIdx = strcmp(projTbl.experimentId,animals{a});
%     ages(a) = unique(projTbl.age(aniIdx));
%     cd(fullfile(dataFold,'Ephys',animals{a}))
%     folders = dir;
%     fileIdx = find(contains({folders.name},'MMM'));
%     if isempty(fileIdx)
%         sumStats = vertcat(sumStats,{[],[],[]});
%         continue
%     end
%     mergeName{a,1} = folders(fileIdx).name;
%     mergeId = mergeName{a}(12:end);
%     load(fullfile(dataFold,'Ephys',animals{a},mergeName{a},[mergeName{a} '_id.mat']))
%     probeId = find(strcmp({id.probes.area},'PSS'));
%     mergeName{a,1} = [mergeName{a,1} '_p' num2str(probeId)];
%     disp(['generating sumStats for ' mergeName{a,1}])
%     [out] = plotMerge(animals{a}, mergeId, probeId, dataFold, 0);
%     sumStats = vertcat(sumStats,out);
% end

load(fullfile(dataFold,'dataSets','cooling',proj,'matchedSU',[proj '_matchedSUdataSet.mat']))

%% Organize Data

preDat = vertcat(sumStats{:,1}); goodId(:,1) = screenUnits(preDat,'SU');
coolDat = vertcat(sumStats{:,2}); goodId(:,2) = screenUnits(coolDat,'SU');
postDat = vertcat(sumStats{:,3}); goodId(:,3) = screenUnits(postDat,'SU');
for u = 1:height(preDat)
    uAge(u,1) = ages(strcmp(animals,preDat.exptName{u}(1:5)));
end

LDR = [];
RP = [];
AGE = [];
for a = 1:length(animals)
    if isempty(sumStats{a,1})
        continue
    end
    dat = sumStats(a,:);
    nU = height(dat{1});

    goodUnitsCntrl = screenUnits(dat{1},anaMode);
    goodUnitsCool = screenUnits(dat{2},anaMode);

    cntrlLDR = dat{1}.ldr; cntrlLDR(~goodUnitsCntrl) = nan;
    coolLDR = dat{2}.ldr; coolLDR(~goodUnitsCool) = nan;
    LDR = vertcat(LDR,[cntrlLDR coolLDR]);

    cntrlrPref = dat{1}.rPref;
    coolrPref = dat{2}.rPref;
    RP = vertcat(RP,[cntrlrPref coolrPref]);

    AGE = vertcat(AGE,repmat(ages(a),nU,1));

    for u = 1:nU
        meanCntrlR = mean(dat{1}.response{u},'omitnan');
        semCntrlR = std(dat{1}.response{u},'omitnan')/sqrt(size(dat{1}.response{u},1));
        cCntrl = dat{1}.condition{u}(strcmp(dat{1}.paramKey{u},'ori'),:);
        meanCoolR = mean(dat{2}.response{u},'omitnan');
        semCoolR = std(dat{2}.response{u},'omitnan')/sqrt(size(dat{2}.response{u},1));
        cCool = dat{2}.condition{u}(strcmp(dat{2}.paramKey{u},'ori'),:);
        meanPostR = mean(dat{3}.response{u},'omitnan');
        semPostR = std(dat{3}.response{u},'omitnan')/sqrt(size(dat{3}.response{u},1));
        cPost = dat{3}.condition{u}(strcmp(dat{3}.paramKey{u},'ori'),:);

        [c,cntrlTC{a}(u,:),alignIdx{a}(u,:,1)] = alignDirTuning(cCntrl,meanCntrlR);
%         [c,coolTC{a}(u,:),alignIdx{a}(u,:,2)] = alignDirTuning(cCool,meanCoolR);
%         [c,postCoolTC{a}(u,:),alignIdx{a}(u,:,3)] = alignDirTuning(cPost,meanPostR);
%         semCoolTC{a}(u,:) = semCoolR(alignIdx{a}(u,:,2));
%         semPostTC{a}(u,:) = semPostR(alignIdx{a}(u,:,3));

        semCntrlTC{a}(u,:) = semCntrlR(alignIdx{a}(u,:,1));
        coolTC{a}(u,:) = meanCoolR(alignIdx{a}(u,:,1));
        semCoolTC{a}(u,:) = semCoolR(alignIdx{a}(u,:,1));
        postCoolTC{a}(u,:) = meanPostR(alignIdx{a}(u,:,1));
        semPostTC{a}(u,:) = semPostR(alignIdx{a}(u,:,1));
    end

end
TC(:,:,1) = vertcat(cntrlTC{:});
TCsem(:,:,1) = vertcat(semCntrlTC{:});
TC(:,:,2) = vertcat(coolTC{:});
TCsem(:,:,2) = vertcat(semCoolTC{:});
TC(:,:,3) = vertcat(postCoolTC{:});
TCsem(:,:,3) = vertcat(semPostTC{:});
TC(TC<0) = 0;
for tId = 1:3
    maxR = max(TC(:,:,1),[],2);
    maxR(maxR==0) = 1;
    TCnorm(:,:,tId) = TC(:,:,tId)./maxR;
end
dTCnorm = (TC(:,:,2)-TC(:,:,1))./(TC(:,:,1)+TC(:,:,2));
% dTCnorm = (TCnorm(:,:,2)-TCnorm(:,:,1))./(TCnorm(:,:,1)+TCnorm(:,:,2));
% dTCnorm = (TC(:,:,2)./(TC(:,:,1)+TC(:,:,2)));


% Compute PSTH for each unit
for j = 1:3
    if j == 1
        data = preDat;
    elseif j == 2
        data = coolDat;
    elseif j == 3
        data = postDat;
    end
    binSize = 0.05;
    edges = -1:binSize:2;
    for u = 1:height(data)
        prefTrials = ismember(data.spkTimes{u}(2,:),  find(data.condition{u}(contains(data.paramKey{u},'ori'),:) == data.oriPref(u)) );
    train(u,:,j) = histcounts(data.spkTimes{u}(1,:),edges)/(max(data.fr(u).trialNum,[],'all')*binSize);
    % train(u,:,j) = train(u,:,j)/max(train(u,:,j));
    end
    %[s,i] = sort(data.latency(~isnan(data.latency)));
    %[s,i] = sort(uAge);
    % train = train(i,:,j);
    % figure; imagesc(train(:,:,j))
end


%% Plot

% figure; hold on
% LDRdiff = LDR(:,2)-LDR(:,1);
% plot(AGE,LDRdiff,'ko','MarkerFaceColor','k')
% uAges = unique(ages);
% for a = 1:length(uAges)
%     x(a) = mean(LDRdiff(AGE==uAges(a)),'omitnan');
%     sem(a) = std(LDRdiff(AGE==uAges(a)),'omitnan')/sqrt(sum(AGE==uAges(a)));
% end
% % plot(uAges,x,'ro')
% errorbar(uAges,x,sem,'r','LineWidth',2,'LineStyle','none')
% yline(0,'k--')
% xlabel('age')
% ylabel('delta Ldir')
% clear x
% box on
% 
% figure;hold on
% ageGroups = {[28 33],[34 40],[41 300]};
% clrs = {[1 0 0], [0.5 0 0], [0 0 0]};
% plot([0 1],[0 1],'k')
% for ag = 1:length(ageGroups)
%     idx = AGE>=ageGroups{ag}(1) & AGE<=ageGroups{ag}(2);
%     plot(LDR(idx,1),LDR(idx,2),'.','Color',clrs{ag},'MarkerSize',10)
% end
% xlabel('control LDR')
% ylabel('cooled LDR')
% 
% 
% 
% figure;hold on
% SI = (RP(:,2) - RP(:,1))./(RP(:,2) + RP(:,1));
% plot(AGE,SI,'ko','MarkerSize',3,'MarkerFaceColor','k')
% % for a = 1:length(uAges)
% %     x(a) = mean(SI(AGE==uAges(a)),'omitnan');
% % end
% % plot(uAges,x,'k*')
% ageBins = min(uAges):3:max(uAges);
% for ab = 1:length(ageBins)-1
%     ageIdx = AGE>ageBins(ab) & AGE<=ageBins(ab+1);
%     binnedAgeSImean(ab) = mean(SI(ageIdx),'omitnan');
%     binnedAgeSIsem(ab) = std(SI(ageIdx),'omitnan')/sqrt(sum(ageIdx));
%     binnedAgeX(ab) = mean(ageBins(ab:ab+1));
% end
% nanIdx = isnan(binnedAgeSImean);
% binnedAgeSImean = binnedAgeSImean(~nanIdx);
% binnedAgeSIsem = binnedAgeSIsem(~nanIdx);
% binnedAgeX = binnedAgeX(~nanIdx);
% plot(binnedAgeX,binnedAgeSImean,'r','LineWidth',2)
% plot(repmat(binnedAgeX,2,1),binnedAgeSImean+([-1;1]*binnedAgeSIsem),'r','LineWidth',2)
% xlabel('age')
% ylabel('SI')
% box on
% 
% ageGroups = {[min(AGE) 32],[33 36],[37 42],[43 max(AGE)]};
% clrs = {[0.4660 0.6740 0.1880],[0.9290 0.6940 0.1250],[0.8500 0.3250 0.0980],[0.6350 0.0780 0.1840]};
% figure;hold on
% for ag = 1:length(ageGroups)
%     cdf = cdfplot(SI(AGE>=ageGroups{ag}(1) & AGE<=ageGroups{ag}(2)));
%     cdf.Color = clrs{ag};
%     cdf.LineWidth = 2;
%     lbl{ag} = ['P' num2str(ageGroups{ag}(1)) '-' num2str(ageGroups{ag}(2))];
% end
% legend(lbl);
% xlabel('SI')
% ylabel('percentile')
% 
% 
% figure; hold on
% plot(LDRdiff,SI,'k.')
% xline(0,'--')
% xlabel('delta Ldir')
% yline(0,'k--')
% ylabel('SI')



clrs = hsv(length(c)-1); clrs = [clrs;clrs(1,:)];

figure;hold on
colororder({'k','r'})
yyaxis left
x = c;
y = mean(TCnorm(:,:,1),'omitnan');
err = std(TCnorm(:,:,1),'omitnan')/sqrt(size(TCnorm,1));
errorbar(x,y,err,'k-','LineWidth',2)
y = mean(TCnorm(:,:,2),'omitnan');
err = std(TCnorm(:,:,2),'omitnan')/sqrt(size(TCnorm,1));
errorbar(x,y,err,'c-','LineWidth',2)
ylabel('response (Hz) norm to rPref cntrl')
scatter(x,0,500,clrs,'.')
yyaxis right
y = mean(dTCnorm,'omitnan');
err = std(dTCnorm,'omitnan')/sqrt(size(dTCnorm,1));
errorbar(x,y,err,'r-','LineWidth',2)
ylabel('SI')
ylim([-1 1])
% y = mean(TCnorm(:,:,2),'omitnan')./mean(TCnorm(:,:,1),'omitnan');
% plot(x,y,'r-','LineWidth',2)
% ylabel('cool/cntrl')
% ylim([0 1])
xticks([-180 -90 0 90 180])
xlabel('dir rel to pref')
box on



figure;hold on
x = mean(TCnorm(:,:,1),'omitnan');
y = mean(TCnorm(:,:,2),'omitnan');
plot(x,y,'k-','LineWidth',2)
scatter(x,y,300,clrs,'.')
plot([0 1],[0 1],'k--')
xlabel('control')
ylabel('V1 cool')
box on



ageGroups = {[29 32],[33 36],[37 40],[41 49],[50 150]};
agClrs = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330]};
figure;hold on
for ag = 1:length(ageGroups)
    idx = uAge>=ageGroups{ag}(1) & uAge<=ageGroups{ag}(2);
    subplot(1,length(ageGroups),ag);hold on
    x = c;
    y = mean(TCnorm(idx,:,1),'omitnan');
    err = std(TCnorm(idx,:,1),'omitnan')/sqrt(sum(idx));
    errorbar(x,y,err,'-','Color',agClrs{ag},'LineWidth',2);
    lbl = ['P' num2str(ageGroups{ag}(1)) '-P' num2str(ageGroups{ag}(2)) '; n=' num2str(sum(idx))];
    y = mean(TCnorm(idx,:,2),'omitnan');
    err = std(TCnorm(idx,:,2),'omitnan')/sqrt(sum(idx));
    errorbar(x,y,err,'--','Color',agClrs{ag},'LineWidth',2)
    xticks([-180 -90 0 90 180])
    title(lbl)
end
box on
clear p lbl



figure; hold on
for ag = 1:length(ageGroups)
    idx = uAge>=ageGroups{ag}(1) & uAge<=ageGroups{ag}(2);
    yline(0,'--')
    x = c;
    y = mean(dTCnorm(idx,:),'omitnan');
    err = std(dTCnorm(idx,:),'omitnan')/sqrt(sum(idx));
    p(ag) = errorbar(x,y,err,'-','LineWidth',2,'Color',agClrs{ag});
%     y = (mean(TCnorm(idx,:,2),'omitnan')-mean(TCnorm(idx,:,1),'omitnan'))./(mean(TCnorm(idx,:,2),'omitnan')+mean(TCnorm(idx,:,1),'omitnan'));
%     p(ag) = plot(x,y,'-','LineWidth',2,'Color',agClrs{ag});
    lbl{ag} = ['P' num2str(ageGroups{ag}(1)) '-P' num2str(ageGroups{ag}(2)) '; n=' num2str(sum(idx))];
    ylabel('SI')
    ylim([-1 1])
    xticks([-180 -90 0 90 180])
    xlabel('dir rel to pref')
end
box on
legend(p,lbl)
clear p lbl



figure; hold on
for ag = 1:length(ageGroups)
    idx = uAge>=ageGroups{ag}(1) & uAge<=ageGroups{ag}(2);
    x = mean(TCnorm(idx,:,1),'omitnan');
    y = mean(TCnorm(idx,:,2),'omitnan');
    p(ag) = plot(x,y,'-','LineWidth',2,'Color',agClrs{ag});
%     scatter(x,y,300,clrs,'.')
    plot([0 1],[0 1],'k--')
    lbl{ag} = ['P' num2str(ageGroups{ag}(1)) '-P' num2str(ageGroups{ag}(2)) '; n=' num2str(sum(idx))];
    xlabel('control')
    ylabel('V1 cool')
end
box on
legend(p,lbl)
clear p lbl



% figure; hold on
% a = TCnorm(:,:,1);
% b = TCnorm(:,:,2);
% c = TCnorm(:,:,3);
% plot(mean(a,'omitnan'))
% plot(mean(b,'omitnan'))
% plot(mean(c,'omitnan'))




% for u = 1:size(TC,1)
% figure('Position',[100 100 600 600]); hold on
% tiledlayout(2,2);
% 
% nexttile; hold on
% errorbar(c,TC(u,:,1),TCsem(u,:,1),'ko-','LineWidth',2)
% errorbar(c,TC(u,:,2),TCsem(u,:,2),'co-','LineWidth',2)
% errorbar(c,TC(u,:,3),TCsem(u,:,3),'bo-','LineWidth',2)
% xticks([-180 -90 0 90 180])
% xlabel('dir rel to pref')
% box on
% 
% nexttile; hold on
% plot(c,TCnorm(u,:,1),'ko-','LineWidth',2)
% plot(c,TCnorm(u,:,2),'co-','LineWidth',2)
% plot(c,TCnorm(u,:,3),'bo-','LineWidth',2)
% xticks([-180 -90 0 90 180])
% xlabel('dir rel to pref')
% box on
% 
% nexttile; hold on
% plot(edges(2:end),train(u,:,1),'k','LineWidth',2)
% plot(edges(2:end),train(u,:,2),'c','LineWidth',2)
% plot(edges(2:end),train(u,:,3),'b','LineWidth',2)
% box on
% 
% nexttile; hold on
% plot(c,dTCnorm(u,:),'ro-','LineWidth',2)
% xticks([-180 -90 0 90 180])
% xlabel('dir rel to pref')
% box on
% 
% end


















































% % % clear all
% % % % close all
% % % 
% % % if ispc
% % % %     dataFold = 'D:\data'; 
% % % %     dataFold = 'C:\Users\brand\Documents\data';
% % % %     dataFold = 'F:\Brandon\data';
% % % dataFold = 'F:\Brandon\VSS2024\data';
% % % elseif ismac
% % % %     dataFold = '/Volumes/Lab drive/Brandon/data';
% % %     dataFold = '/Volumes/Lab drive/Brandon/VSS2024/data';
% % % %     dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
% % % end
% % % physDir = fullfile(dataFold,'Ephys');
% % % dsDir = fullfile(dataFold,'dataSets');
% % % 
% % % anaMode = 'SU';
% % % for i = 1:3
% % %     figure; hold on
% % %     if i == 1
% % %         load(fullfile(dataFold,'dataSets',['P33to34Merge' anaMode 'Dat.mat']))
% % %         clr = [0.9290 0.6940 0.1250]; %yellow
% % %         ttl = 'P33-34';
% % %     elseif i == 2
% % %         load(fullfile(dataFold,'dataSets',['P36to37Merge' anaMode 'Dat.mat']))
% % %         clr = [0.8500 0.3250 0.0980]; %orange
% % %         ttl = 'P36-37';
% % %     elseif i == 3
% % %         load(fullfile(dataFold,'dataSets',['adultMerge' anaMode 'Dat.mat']))
% % %         clr = [0.6350 0.0780 0.1840]; %red
% % %         ttl = 'adult';
% % %     end
% % %     fileID = vertcat(dat.tbl.fileID{:});
% % %     condID = fileID(:,2);
% % %     idxA = dat.tbl.goodUnit == 1 & condID == 1;
% % %     idxB = dat.tbl.goodUnit == 1 & condID == 2;
% % % 
% % %     A = dat.tbl(idxA,:).DSI;
% % %     B = dat.tbl(idxB,:).DSI;
% % % 
% % %     A(A>1) = 1;
% % %     B(B>1) = 1;
% % % 
% % %     cdfA = cdfplot(A);
% % %     cdfA.Color = 'k';
% % %     
% % %     cdfB = cdfplot(B);
% % %     cdfB.Color = 'c';
% % % 
% % %     xlim([0 1])
% % %     title(ttl)
% % %     clear dat
% % % end