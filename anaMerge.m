%anaMerge
clear all
close all

proj = 'V1cool_ori';
% dataFold = '/Volumes/NielsenHome2/Brandon/data';
dataFold = 'Y:\Brandon\data';
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
        [x,cntrlTC{a}(u,:),i{a}(u,:)] = alignDirTuning(dat{1}.condition{u}(strcmp(dat{1}.paramKey{u},'ori'),:),meanCntrlR);
        coolTC{a}(u,:) = mean(dat{2}.response{u}(:,i{a}(u,:)),'omitnan');
        postCoolTC{a}(u,:) = mean(dat{3}.response{u}(:,i{a}(u,:)),'omitnan');
    end

end
TC(:,:,1) = vertcat(cntrlTC{:});
TC(:,:,2) = vertcat(coolTC{:});
TC(:,:,3) = vertcat(postCoolTC{:});
TC(TC<0) = 0;
for tId = 1:3
    maxR = max(TC(:,:,tId),[],2);
    maxR(maxR==0) = 1;
    TC(:,:,tId) = TC(:,:,tId)./maxR;
end

%% Plot

figure; hold on
LDRdiff = LDR(:,2)-LDR(:,1);
plot(AGE,LDRdiff,'k.')
uAges = unique(ages);
for a = 1:length(uAges)
    x(a) = mean(LDRdiff(AGE==uAges(a)),'omitnan');
end
plot(uAges,x,'ko')
yline(0,'k--')
xlabel('age')
ylabel('delta Ldir')
clear x

figure;hold on
ageGroups = {[28 33],[34 40],[41 300]};
clrs = {[1 0 0], [0.5 0 0], [0 0 0]};
plot([0 1],[0 1],'k')
for ag = 1:length(ageGroups)
    idx = AGE>=ageGroups{ag}(1) & AGE<=ageGroups{ag}(2);
    plot(LDR(idx,1),LDR(idx,2),'.','Color',clrs{ag},'MarkerSize',10)
end
xlabel('control LDR')
ylabel('cooled LDR')



figure;hold on
SI = (RP(:,2) - RP(:,1))./(RP(:,2) + RP(:,1));
plot(AGE,SI,'ko','MarkerSize',3,'MarkerFaceColor','k')
% for a = 1:length(uAges)
%     x(a) = mean(SI(AGE==uAges(a)),'omitnan');
% end
% plot(uAges,x,'k*')
ageBins = min(uAges):3:max(uAges);
for ab = 1:length(ageBins)-1
    ageIdx = AGE>ageBins(ab) & AGE<=ageBins(ab+1);
    binnedAgeSImean(ab) = mean(SI(ageIdx),'omitnan');
    binnedAgeSIsem(ab) = std(SI(ageIdx),'omitnan')/sqrt(sum(ageIdx));
    binnedAgeX(ab) = mean(ageBins(ab:ab+1));
end
nanIdx = isnan(binnedAgeSImean);
binnedAgeSImean = binnedAgeSImean(~nanIdx);
binnedAgeSIsem = binnedAgeSIsem(~nanIdx);
binnedAgeX = binnedAgeX(~nanIdx);
plot(binnedAgeX,binnedAgeSImean,'r','LineWidth',2)
plot(repmat(binnedAgeX,2,1),binnedAgeSImean+([-1;1]*binnedAgeSIsem),'r','LineWidth',2)
xlabel('age')
ylabel('SI')
box on

ageGroups = {[min(AGE) 32],[33 36],[37 42],[43 max(AGE)]};
clrs = {[0.4660 0.6740 0.1880],[0.9290 0.6940 0.1250],[0.8500 0.3250 0.0980],[0.6350 0.0780 0.1840]};
figure;hold on
for ag = 1:length(ageGroups)
    cdf = cdfplot(SI(AGE>=ageGroups{ag}(1) & AGE<=ageGroups{ag}(2)));
    cdf.Color = clrs{ag};
    cdf.LineWidth = 2;
    lbl{ag} = ['P' num2str(ageGroups{ag}(1)) '-' num2str(ageGroups{ag}(2))];
end
legend(lbl);
xlabel('SI')
ylabel('percentile')


figure; hold on
plot(LDRdiff,SI,'k.')
xline(0,'--')
xlabel('delta Ldir')
yline(0,'k--')
ylabel('SI')


figure;
u = 11;
subplot(1,2,1);hold on
plot(TC(u,:,1),'k')
plot(TC(u,:,2),'c')
subplot(1,2,2);hold on
plot(TC(:,:,1),TC(:,:,2),'.')
plot([0 1],[0 1],'k-')


figure; hold on
a = TC(:,:,1);
b = TC(:,:,2);
c = TC(:,:,3);
plot(mean(a,'omitnan'))
plot(mean(b,'omitnan'))
plot(mean(c,'omitnan'))






































































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