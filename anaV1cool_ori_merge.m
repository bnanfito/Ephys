%anaV1cool_ori_merge
clear all
% close all

proj = 'V1cool_ori';
% dataFold = '/Volumes/NielsenHome2/Brandon/data';
dataFold = 'Y:\Brandon\data';
% dataFold = 'C:\Users\brand\Documents\data';
% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
anaMode = 'SU';
ageGroups = {[28 32],[33 40],[41 80],[81 120]};

%% Build Project Table

% projTbl = getProjectFiles(proj,1,'age','recSite','penNr','priorMFlag','priorDescr',...
%                                        'duringMFlag','manipDescr','manipDetail',...
%                                        'looperNameCond1','looperNameCond2',...
%                                        'looperNameCond3','looperNameCond4',...
%                                        'looperNameCond5');
% animals = unique(projTbl.experimentId);
% sumStats = [];
% for a = 1:length(animals)
% %     if ~strcmp(animals{a},'febn2')
% %         continue
% %     end
%     aniIdx = strcmp(projTbl.experimentId,animals{a});
%     ages(a) = unique(projTbl.age(aniIdx));
%     cd(fullfile(dataFold,'Ephys',animals{a}))
%     folders = dir;
%     fileIdx = find(contains({folders.name},'MMM'));
%     if isempty(fileIdx)
%         sumStats = vertcat(sumStats,{[],[],[]});
%         continue
%     end
% 
%     out = cell(1,3);
%     for fId = 1:length(fileIdx)
%         mergeName{a,fId} = folders(fileIdx(fId)).name;
%         mergeId = mergeName{a,fId}(12:end);
%         load(fullfile(dataFold,'Ephys',animals{a},mergeName{a,fId},[mergeName{a,fId} '_id.mat']))
%         probeId = find(strcmp({id.probes.area},'PSS'));
%         disp(['generating sumStats for ' mergeName{a,fId}])
%     
% %         splitIntan(fullfile(dataFold,'Ephys'),animals{a},mergeId,probeId,'BRN')
%         %Load merge info
%         load(fullfile(dataFold,'Ephys',animals{a},mergeName{a,fId},[mergeName{a,fId} '_mergeInfo.mat']))
%         nFiles = length(mergeInfo.files);
%         for f = 1:nFiles
%             exptName{f,1} = [animals{a} '_' mergeInfo.files{f}];
%             load(fullfile(dataFold,'Ephys',animals{a},exptName{f},[exptName{f} '_trialInfo.mat']))
%             dom = trialInfo.dom;
%             isOri = (length(dom)==1 & sum(strcmp(dom,'ori'))==1) |...
%                     (length(dom)==2 & sum(strcmp(dom,'ori'))==1 & sum(contains(dom,'size'))==1 );
%             if isOri
%                 out_tmp = anaOri(animals{a},exptName{f}(8:10),exptName{f}(12:14),probeId,anaMode,dataFold,0,0,f);
%                 out{f} = vertcat(out{f},out_tmp);
%             end
%         end
%     end
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
for ag = 1:length(ageGroups)
    uAG(uAge>=ageGroups{ag}(1) & uAge<=ageGroups{ag}(2)) = ag;
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
        
%         if e == 1
            [c_aligned,tc_aligned(u,:,e),alignIdx(u,:)] = alignDirTuning(conds,tc_temp);
%         else
%             tc_aligned = tc_temp(alignIdx(u,:));
%         end
        tc_sem_algined = tc_sem_temp(alignIdx(u,:));

        maxR = max(tc_aligned(u,:,e));
        tc_norm(u,:,e) = tc_aligned(u,:,e)/maxR;

        tc(u,:,e) = tc_temp;
        tc_sem(u,:,e) = tc_sem_temp;

        c = unique(abs(c_aligned));
        for i = 1:length(c)
            tc_alignedH(u,i,e) = mean(tc_aligned(u,abs(c_aligned)==c(i),e),'omitnan');
            tc_normH(u,i,e) = mean(tc_norm(u,abs(c_aligned)==c(i),e),'omitnan');
        end

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

%% data table

r1 = tc_alignedH(:,:,1);
r1_norm = tc_normH(:,:,1);
r2 = tc_alignedH(:,:,2);
r2_norm = tc_normH(:,:,2);
c = repmat(c,height(dat{1}),1);

d = table();
d.r1 = r1(:);
d.r1L = log10(d.r1+1);
d.r1_norm = r1_norm(:);
d.r2 = r2(:);
d.r2L = log10(d.r2+1);
d.r2_norm = r2_norm(:);
d.si = (d.r2-d.r1)./(d.r2+d.r1);
d.siL = (d.r2L-d.r1L)./(d.r2L+d.r1L);
d.logRatio = log(d.r2./d.r1);
d.ori = c(:);

%% Plot

figure('Position',[100 100 1000 1000]);
mSize = 7;

count = 0;
exUs = [4 25];
exUclrs = {[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880]};
agShapes = {'+','square','diamond','^'};
for u = exUs
    subplot(3,4,1+(count*2)); hold on
    errorbar(conds,tc(u,:,1),tc_sem(u,:,1),'k','LineWidth',1)
    errorbar(conds,tc(u,:,2),tc_sem(u,:,2),'c','LineWidth',1)
    % errorbar(conds,tc(u,:,3),tc_sem(u,:,3),'b')
    plot(conds,tc(u,:,2)*(rPref(u,1)/rPref(u,2)),'c--','LineWidth',1)
    % legend('cntrl','cool','scaled cool')
    xticks([0 90 180 270])
    xlim([0 360])
    xlabel('dir of motion (deg)')
    ylabel('firing rate')
    box on
    axis square
    
    subplot(3,4,2+(count*2)); hold on
    countTrials = 0;
    for i = 1:3
        if i==1
            clr = 'k';
        elseif i==2
            clr = 'c';
        elseif i==3
            clr = 'k';
        end
        x = dat{i}.spkTimes{u}(1,:);
        y = dat{i}.spkTimes{u}(2,:);
        y = y+countTrials;
        plot(x,y,[clr '.'],'MarkerSize',4)
        countTrials = countTrials+max(dat{i}.fr(u).trialNum,[],'all');
    end
    xlim([-1 2])
    xlabel('time (sec)')
    ylim([0 countTrials+1])
    ylabel('trial number')
    patch([0 1 1 0],[0 0 countTrials+1 countTrials+1],'k','EdgeColor','none','FaceAlpha',0.2)
    box on
    axis square
    
    count = count+1;
end

subplot(3,3,4); hold on
x = ldr(:,2);
y = ldr(:,1);
% plot(x,y,'k.','MarkerSize',4)
for ag = unique(uAG)
    p(ag) = plot(x(uAG==ag),y(uAG==ag),['k' agShapes{ag}],'MarkerSize',mSize,'LineWidth',1);
    lbl{ag} = ['P' num2str(ageGroups{ag}(1)) '-P' num2str(ageGroups{ag}(2)) '; n=' num2str(sum(uAG==ag))];
end
plot(x(exUs(1)),y(exUs(1)),'o','MarkerSize',mSize,'Color',exUclrs{1},'LineWidth',1)
plot(x(exUs(2)),y(exUs(2)),'o','MarkerSize',mSize,'Color',exUclrs{2},'LineWidth',1)
plot([0 1],[0 1],'k--')
xlabel('V1 cooled')
ylabel('control')
title('LDR')
legend(p,lbl,'Location','northwest')
box on
axis square

subplot(3,3,5); hold on
x = dsi(:,2);
y = dsi(:,1);
% plot(x,y,'k.','MarkerSize',4)
for ag = 1:length(ageGroups)
    plot(x(uAG==ag),y(uAG==ag),['k' agShapes{ag}],'MarkerSize',mSize,'LineWidth',1)
end
plot(x(exUs(1)),y(exUs(1)),'o','MarkerSize',mSize,'Color',exUclrs{1},'LineWidth',1)
plot(x(exUs(2)),y(exUs(2)),'o','MarkerSize',mSize,'Color',exUclrs{2},'LineWidth',1)
plot([0 1],[0 1],'k--')
xlabel('V1 cooled')
ylabel('control')
title('DSI')
sgtitle('tuning metrics')
box on
axis square

subplot(3,3,6); hold on
x = bwS(:,2);
y = bwS(:,1);
lims = [10 60];
out = x>lims(2)|y>lims(2);
% plot(x,y,'k.','MarkerSize',4)
for ag = 1:length(ageGroups)
    plot(x(uAG==ag),y(uAG==ag),['k' agShapes{ag}],'MarkerSize',mSize,'LineWidth',1)
end
x(x>lims(2)) = lims(2);
y(y>lims(2)) = lims(2);
for ag = 1:length(ageGroups)
    plot(x(uAG==ag & out'),y(uAG==ag & out'),['r' agShapes{ag}],'MarkerSize',mSize,'LineWidth',1)
end
plot(x(exUs(1)),y(exUs(1)),'o','MarkerSize',mSize,'Color',exUclrs{1},'LineWidth',1)
plot(x(exUs(2)),y(exUs(2)),'o','MarkerSize',mSize,'Color',exUclrs{2},'LineWidth',1)
plot([0 200],[0 200],'k--')
xlabel('V1 cooled') 
xlim(lims)
ylim(lims)
ylabel('control')
title('bandwidth (smooth)')
box on
axis square


subplot(3,2,5); hold on
colororder({'k','c'})
x = table();
x.ori = repmat(d.ori,2,1);
x.manip = vertcat(repmat({'control'},height(d),1),repmat({'cool'},height(d),1));
% x.r = vertcat(d.r1,d.r2)+1;
x.r = vertcat(d.r1L,d.r2L);
boxchart(categorical(x.ori),x.r,'GroupByColor',x.manip,'Notch','on')
% ylabel('R')
ylabel('log(R+1)')
% set(gca,'YScale','log')
xlabel('Angular disparity (relative to pref)')
box on

subplot(3,2,6); hold on
boxchart(categorical(d.ori),d.si,'notch','on')
ylabel('SI')
% boxchart(categorical(d.ori),d.siL,'notch','on')
% ylabel('SI (log)')
yline(0,'k--')
xlabel('Angular disparity (relative to pref)')
box on





% close all
% clear dirTuning
% e = 1;
% for u = 1:height(dat{e})
% fr = dat{e}.response{u}; fr(fr<0) = 0;
% dirTuning(u) = getDirTuning(fr,conds,1);
% sgtitle([dat{e}.exptName{u} ' ' dat{e}.uInfo{u} num2str(dat{e}.uID(u))])
% saveas(gcf,[dat{e}.uInfo{u} num2str(u)])
% end
% save('dirTuning','dirTuning')



% u = 25;
% figure; tiledlayout('flow')
% nexttile; hold on
% errorbar(conds,tc(u,:,1),tc_sem(u,:,1),'k','LineWidth',2)
% errorbar(conds,tc(u,:,2),tc_sem(u,:,2),'c','LineWidth',2)
% % errorbar(conds,tc(u,:,3),tc_sem(u,:,3),'b')
% plot(conds,tc(u,:,2)*(rPref(u,1)/rPref(u,2)),'c--','LineWidth',2)
% % legend('cntrl','cool','scaled cool')
% xticks([0 90 180 270])
% xlim([0 360])
% xlabel('dir of motion (deg)')
% ylabel('firing rate')
% box on
% axis square
% nexttile; hold on
% countTrials = 0;
% for i = 1:3
%     if i==1
%         clr = 'k';
%     elseif i==2
%         clr = 'c';
%     elseif i==3
%         clr = 'k';
%     end
%     x = dat{i}.spkTimes{u}(1,:);
%     y = dat{i}.spkTimes{u}(2,:);
%     y = y+countTrials;
%     scatter(x,y,[clr '.'])
%     countTrials = countTrials+max(dat{i}.fr(u).trialNum,[],'all');
% end
% xlim([-1 2])
% ylim([0 countTrials+1])
% patch([0 1 1 0],[0 0 countTrials+1 countTrials+1],'k','EdgeColor','none','FaceAlpha',0.2)
% box on
% axis square
% sgtitle(['example unit: ' dat{1}.exptName{u} ' ' dat{1}.uInfo{u} num2str(dat{1}.uID(u))])

% nexttile; hold on
% plot(c_aligned,tc_norm(u,:,1),'k')
% plot(c_aligned,tc_norm(u,:,2),'c')
% %     plot(c_aligned,tc_norm(u,:,3),'b')
% xticks([-180 -90 0 90 180])
% xlabel('dir rel to pref')
% ylabel('normalized response')
% title('normalized example unit')
% box on
% axis square
% nexttile; hold on
% y = mean(tc_norm(:,:,1),'omitnan');
% sem = std(tc_norm(:,:,1),'omitnan')/sqrt(size(tc_norm,1));
% errorbar(c_aligned,y,sem,'k')
% y = mean(tc_norm(:,:,2),'omitnan');
% sem = std(tc_norm(:,:,2),'omitnan')/sqrt(size(tc_norm,1));
% errorbar(c_aligned,y,sem,'c')
% xticks([-180 -90 0 90 180])
% xlabel('dir rel to pref')
% ylabel('normalized response')
% title('normalized population')
% box on
% axis square
% nexttile; hold on
% plot(bins(2:end),psth(u,:,1),'k')
% plot(bins(2:end),psth(u,:,2),'c')
% %     plot(bins(2:end),psth(u,:,3),'b')
% xlabel('time (sec)')
% ylabel('firing rate')
% box on
% axis square
% nexttile; hold on
% plot(tc_norm(u,:,1),tc_norm(u,:,2),'ko-')
% plot([0 1],[0 1],'k--')
% xlabel('control norm response')
% ylabel('v1 cooled norm response')
% box on
% axis square
% nexttile;hold on
% plot(mean(tc_norm(:,:,1),'omitnan'),mean(tc_norm(:,:,2),'omitnan'),'k-o')
% plot([0 1],[0 1],'k--')
% xlabel('control response')
% ylabel('V1 cooled response')
% box on
% axis square

% figure;tiledlayout(2,2)
% nexttile;hold on
% x = bw(:,2);
% y = bw(:,1);
% scatter(x,y,'k.')
% plot([0 max(bw,[],'all')],[0 max(bw,[],'all')],'k--')
% xlabel('cool')
% ylabel('cntrl')
% title('bandwidth')
% box on
% axis square
% nexttile;hold on
% x = bwS(:,2);
% y = bwS(:,1);
% scatter(x,y,'k.')
% plot([0 max(bwS,[],'all')],[0 max(bwS,[],'all')],'k--')
% xlabel('cool')
% ylabel('cntrl')
% title('bandwidth (smooth)')
% box on
% axis square
% nexttile;hold on
% x = ldr(:,2);
% y = ldr(:,1);
% scatter(x,y,'k.')
% plot([0 1],[0 1],'k--')
% xlabel('cool')
% ylabel('cntrl')
% title('LDR')
% box on
% axis square
% nexttile;hold on
% x = dsi(:,2);
% y = dsi(:,1);
% scatter(x,y,'k.')
% plot([0 1],[0 1],'k--')
% xlabel('cool')
% ylabel('cntrl')
% title('DSI')
% sgtitle('tuning metrics')
% box on
% axis square

% figure;tiledlayout(2,2)
% nexttile;hold on
% x = rPref(:,1);
% y = rPref(:,2);
% scatter(x,y,'k.')
% plot([0 max([x y],[],'all')],[0 max([x y],[],'all')],'k--')
% title('rPref')
% xlabel('cntrl')
% ylabel('cool')
% box on
% axis square
% nexttile;hold on
% x = rNull(:,1);
% y = rNull(:,2);
% scatter(x,y,'k.')
% plot([0 max([x y],[],'all')],[0 max([x y],[],'all')],'k--')
% title('rNull')
% xlabel('cntrl')
% ylabel('cool')
% box on
% axis square
% nexttile;hold on
% cdf = cdfplot(rPref(:,1));
% cdf.Color = 'k';
% cdf = cdfplot(rPref(:,2));
% cdf.Color = 'c';
% xlabel('rPref')
% ylabel('proportion')
% yticks([0 0.5 1])
% yline(0.5,'k--')
% grid off
% box on
% axis square
% nexttile;hold on
% cdf = cdfplot(rNull(:,1));
% cdf.Color = 'k';
% cdf = cdfplot(rNull(:,2));
% cdf.Color = 'c';
% xlabel('rNull')
% ylabel('proportion')
% yticks([0 0.5 1])
% yline(0.5,'k--')
% grid off
% box on
% axis square

% figure; tiledlayout(2,2)
% nexttile; hold on
% x = rNull(:,2);
% y = scldNull;
% dRnull_scl = x-y;
% scatter(x,y,'k.')
% plot([0 max([x y],[],'all')],[0 max([x y],[],'all')],'k--')
% title('rNull')
% xlabel('cool')
% ylabel('scaled')
% box on
% axis square
% nexttile; hold on
% idx = c_aligned==180;
% x = tc_norm(:,idx,2);
% y = tc_norm(:,idx,1);
% dRnull_norm = x-y;
% scatter(x,y,'k.')
% plot([0 1],[0 1],'k--')
% title('normalized null')
% xlabel('cool')
% ylabel('cntrl')
% box on
% axis square
% nexttile; hold on
% histogram(dRnull_scl)
% [p,h] = signrank(dRnull_scl);
% xline(mean(dRnull_scl,'omitnan'),'r')
% xline(median(dRnull_scl,'omitnan'),'g')
% ylabel('count')
% xlabel('cool null - scaled null')
% title(['signrank: p=' num2str(p)])
% box on
% axis square
% nexttile; hold on
% histogram(dRnull_norm)
% [p,h] = signrank(dRnull_norm);
% xline(mean(dRnull_norm,'omitnan'),'r')
% xline(median(dRnull_norm,'omitnan'),'g')
% ylabel('count')
% xlabel('delta norm. rNull (cool-cntrl)')
% title(['signrank: p=' num2str(p)])
% box on
% axis square

% figure;hold on
% x = scl;
% y = scl_null;
% sc = scatter(x,y,'k.');
% dt = dataTipTextRow('uID',dat{1}.uID);
% sc.DataTipTemplate.DataTipRows(end+1) = dt;
% plot([0 2],[0 2],'k--')
% xline(1,'k--')
% yline(1,'k--')
% xlabel('scale pref')
% ylabel('scale null')
