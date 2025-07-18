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

% compile sumStat tables by phase of experiment (1:cntrl, 2:cool, 3:post)
dat{1} = vertcat(sumStats{:,1});
dat{2} = vertcat(sumStats{:,2});
dat{3} = vertcat(sumStats{:,3});

% determine which units pass inclusion criteria for each phase of the expt
goodId(:,1) = screenUnits(dat{1},'SU');
goodId(:,2) = screenUnits(dat{2},'SU');
goodId(:,3) = screenUnits(dat{3},'SU');

keepIdx = goodId(:,1)&goodId(:,2);
dat{1} = dat{1}(keepIdx,:);
dat{2} = dat{2}(keepIdx,:);
dat{3} = dat{3}(keepIdx,:);

% for each phase, extract basic information like age, tuning metrics, etc
% and store that info in vectors
nU = height(dat{1});
uAge = nan(nU,1);
ldr = nan(nU,3);
dsi = nan(nU,3);
bw = nan(nU,3);
bwS = nan(nU,3);
late = nan(nU,3);
rPref = nan(nU,3);
rNull = nan(nU,3);
for u = 1:nU
    uAge(u) = ages(strcmp(animals,dat{1}.exptName{u}(1:5)));
end
for e = 1:3
%     ldr(:,e) = dat{e}.ldr;
    late(:,e) = dat{e}.latency;
%     rPref(:,e) = dat{e}.rPref;
    for u = 1:nU

        rMat = dat{e}.response{u};
        rMat(rMat<0) = 0;
        conds = dat{e}.condition{u}(strcmp(dat{e}.paramKey{u},'ori'),:);

        dirDat = getDirTuning(rMat,conds,0);
        rPref(u,e) = dirDat.prefResp;
        rNull(u,e) = dirDat.nullResp;
        ldr(u,e) = dirDat.Ldir;
        dsi(u,e) = dirDat.DSI;
        bw(u,e) = dirDat.BW;
        bwS(u,e) = dirDat.BWSmooth;

        %tuning curve
        tc_temp = mean(rMat,1,'omitnan');
        tc_sem_temp = std(rMat,[],1,'omitnan')/sqrt(size(rMat,1));
        
        if e == 1
            [c_aligned,tc_aligned,alignIdx(u,:)] = alignDirTuning(conds,tc_temp);
        else
            tc_aligned = tc_temp(alignIdx(u,:));
        end
        tc_sem_algined = tc_sem_temp(alignIdx(u,:));

        maxR = max(tc_aligned);
        tc_norm(u,:,e) = tc_aligned/maxR;

        tc(u,:,e) = tc_temp;
        tc_sem(u,:,e) = tc_sem_temp;

        %PSTH
        binSize = 0.02;
        bins = -1:binSize:2;
        spkTimes = dat{e}.spkTimes{u}(1,:);
        nTrials = max(dat{e}.fr(u).trialNum,[],'all');
        psth(u,:,e) = histcounts(spkTimes,bins)/(nTrials*binSize);

    end
end

%Are pref and null scaled the same during cooling?
scl = rPref(:,2)./rPref(:,1);
scldNull = rNull(:,1).*scl;
scldTc = tc(:,:,1).*scl;
scl_null = rNull(:,2)./rNull(:,1);

%% Plot

u = 4;
figure; tiledlayout(2,3)
nexttile; hold on
errorbar(conds,tc(u,:,1),tc_sem(u,:,1),'k')
errorbar(conds,tc(u,:,2),tc_sem(u,:,2),'c')
%     errorbar(conds,tc(u,:,3),tc_sem(u,:,3),'b')
plot(conds,scldTc(u,:),'c--')
legend('cntrl','cool','scaled cool')
xlabel('dir of motion (deg)')
ylabel('firing rate')
title('example unit')
box on
axis square
nexttile; hold on
plot(c_aligned,tc_norm(u,:,1),'k')
plot(c_aligned,tc_norm(u,:,2),'c')
%     plot(c_aligned,tc_norm(u,:,3),'b')
xlabel('dir rel to pref')
ylabel('normalized response')
title('normalized example unit')
box on
axis square
nexttile; hold on
y = mean(tc_norm(:,:,1),'omitnan');
sem = std(tc_norm(:,:,1),'omitnan')/sqrt(size(tc_norm,1));
errorbar(c_aligned,y,sem,'k')
y = mean(tc_norm(:,:,2),'omitnan');
sem = std(tc_norm(:,:,2),'omitnan')/sqrt(size(tc_norm,1));
errorbar(c_aligned,y,sem,'c')
xticks([-180 -90 0 90 180])
xlabel('dir rel to pref')
ylabel('normalized response')
title('normalized population')
box on
axis square
nexttile; hold on
plot(bins(2:end),psth(u,:,1),'k')
plot(bins(2:end),psth(u,:,2),'c')
%     plot(bins(2:end),psth(u,:,3),'b')
xlabel('time (sec)')
ylabel('firing rate')
box on
axis square
nexttile; hold on
plot(tc_norm(u,:,1),tc_norm(u,:,2),'ko-')
plot([0 1],[0 1],'k--')
xlabel('control norm response')
ylabel('v1 cooled norm response')
box on
axis square
nexttile;hold on
plot(mean(tc_norm(:,:,1),'omitnan'),mean(tc_norm(:,:,2),'omitnan'),'k-o')
plot([0 1],[0 1],'k--')
xlabel('control response')
ylabel('V1 cooled response')
box on
axis square

figure;tiledlayout(2,2)
nexttile;hold on
scatter(bw(:,1),bw(:,2),'k.')
plot([0 max(bw,[],'all')],[0 max(bw,[],'all')],'k--')
xlabel('control bandwidth')
ylabel('v1 cooled bandwidth')
box on
axis square
nexttile;hold on
scatter(bwS(:,1),bwS(:,2),'k.')
plot([0 max(bwS,[],'all')],[0 max(bwS,[],'all')],'k--')
xlabel('control bandwidth (smooth)')
ylabel('v1 cooled bandwidth (smooth)')
box on
axis square
nexttile;hold on
scatter(ldr(:,1),ldr(:,2),'k.')
plot([0 1],[0 1],'k--')
xlabel('control ldr')
ylabel('v1 cooled ldr')
box on
axis square
nexttile;hold on
scatter(dsi(:,1),dsi(:,2),'k.')
plot([0 1],[0 1],'k--')
xlabel('control dsi')
ylabel('v1 cooled dsi')
sgtitle('tuning metrics')
box on
axis square

figure;tiledlayout(2,2)
nexttile;hold on
x = rPref(:,1);
y = rPref(:,2);
scatter(x,y,'k.')
plot([0 max([x y],[],'all')],[0 max([x y],[],'all')],'k--')
title('rPref')
xlabel('cntrl')
ylabel('cool')
box on
axis square
nexttile;hold on
x = rNull(:,1);
y = rNull(:,2);
scatter(x,y,'k.')
plot([0 max([x y],[],'all')],[0 max([x y],[],'all')],'k--')
title('rNull')
xlabel('cntrl')
ylabel('cool')
box on
axis square
nexttile;hold on
cdf = cdfplot(rPref(:,1));
cdf.Color = 'k';
cdf = cdfplot(rPref(:,2));
cdf.Color = 'c';
xlabel('rPref')
ylabel('proportion')
yticks([0 0.5 1])
yline(0.5,'k--')
grid off
box on
axis square
nexttile;hold on
cdf = cdfplot(rNull(:,1));
cdf.Color = 'k';
cdf = cdfplot(rNull(:,2));
cdf.Color = 'c';
xlabel('rNull')
ylabel('proportion')
yticks([0 0.5 1])
yline(0.5,'k--')
grid off
box on
axis square

figure; tiledlayout(2,2)
nexttile; hold on
x = scldNull;
y = rNull(:,2);
dRnull_scl = y-x;
scatter(x,y,'k.')
plot([0 max([x y],[],'all')],[0 max([x y],[],'all')],'k--')
title('rNull')
xlabel('scaled')
ylabel('cool')
box on
axis square
nexttile; hold on
idx = c_aligned==180;
x = tc_norm(:,idx,1);
y = tc_norm(:,idx,2);
dRnull_norm = y-x;
scatter(x,y,'k.')
plot([0 1],[0 1],'k--')
title('normalized null')
xlabel('cntrl')
ylabel('cool')
box on
axis square
nexttile; hold on
histogram(dRnull_scl)
[p,h] = signrank(dRnull_scl);
xline(mean(dRnull_scl,'omitnan'),'r')
xline(median(dRnull_scl,'omitnan'),'g')
ylabel('count')
xlabel('cool null - scaled null')
title(['signrank: p=' num2str(p)])
box on
axis square
nexttile; hold on
histogram(dRnull_norm)
[p,h] = signrank(dRnull_norm);
xline(mean(dRnull_norm,'omitnan'),'r')
xline(median(dRnull_norm,'omitnan'),'g')
ylabel('count')
xlabel('delta norm. rNull (cool-cntrl)')
title(['signrank: p=' num2str(p)])
box on
axis square

figure;hold on
x = scl;
y = scl_null;
sc = scatter(x,y,'k.');
dt = dataTipTextRow('uID',dat{1}.uID);
sc.DataTipTemplate.DataTipRows(end+1) = dt;
plot([0 2],[0 2],'k--')
xline(1,'k--')
yline(1,'k--')
xlabel('scale pref')
ylabel('scale null')
