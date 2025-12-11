%anaV1cool_ori_merge
clear all
close all

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
% dsi = nan(nU,3);
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

        %response matrix
        rMat = dat{e}.response{u};
        rMat(rMat<0) = 0;
        conds = dat{e}.condition{u}(strcmp(dat{e}.paramKey{u},'ori'),:);
        
        %tuning curve
        rMean = mean(rMat,1,'omitnan');
        rSem = std(rMat,[],1,'omitnan')/sqrt(size(rMat,1));

        wHan=hanning(3);
        wHan=wHan/sum(wHan);
        rWrap = [rMean rMean rMean];
        rSmooth = conv(rWrap,wHan,'same');
        rSmooth = rSmooth((1:length(rMean))+length(rMean));

        [dsi(u,e),~,~,DSI(u,e,:),R_PREF(u,e,:),C_PREF(u,e,:),dsi_s(u,e),~,~,DSI_sRep(u,e,:)] = compute_dsi(conds,rMat,0);
        mV = meanvec(conds,rMean);
        ldr(u,e) = mV.ldr;
        angDir(u,e) = mV.angDir;
        mV_s = meanvec(conds,rSmooth);
        ldr_s(u,e) = mV_s.ldr;
        angDir_s(u,e) = mV_s.angDir;
        nReps = 5;
        for rep = 1:nReps
            rRep = rMat(rep,:);

            wHan=hanning(3);
            wHan=wHan/sum(wHan);
            rWrap_rep = [rRep rRep rRep];
            rSmooth_rep = conv(rWrap_rep,wHan,'same');
            rSmooth_rep = rSmooth_rep((1:length(rRep))+length(rRep));

            mV = meanvec(conds,rRep);
            LDR(u,e,rep) = mV.ldr;
            ANG_DIR(u,e,rep) = mV.angDir;
            LOR(u,e,rep) = mV.lor;
            ANG_ORI(u,e,rep) = mV.angOri;

            mV_s = meanvec(conds,rSmooth_rep);
            LDR_s(u,e,rep) = mV_s.ldr;
            ANG_DIR_s(u,e,rep) = mV_s.angDir;
            LOR_s(u,e,rep) = mV_s.lor;
            ANG_ORI_s(u,e,rep) = mV_s.angOri;
        end

        dirDat = getDirTuning(rMat,conds,0);
        rPref(u,e) = dirDat.prefResp;
        rNull(u,e) = dirDat.nullResp;
%         ldr(u,e) = dirDat.Ldir;
%         dsi(u,e) = dirDat.DSI;
        bw(u,e) = dirDat.BW;
        bwS(u,e) = dirDat.BWSmooth;
        
%         if e == 1
            [c_aligned,tc_aligned(u,:,e),alignIdx(u,:)] = alignDirTuning(conds,rMean);
%         else
%             tc_aligned = tc_temp(alignIdx(u,:));
%         end
        tc_sem_algined = rSem(alignIdx(u,:));

        maxR = max(tc_aligned(u,:,e));
        tc_norm(u,:,e) = tc_aligned(u,:,e)/maxR;

        tc(u,:,e) = rMean;
        tc_sem(u,:,e) = rSem;

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
    ciLDR(:,:,e) = confInt(squeeze(LDR(:,e,:))');
end

mLDR = mean(LDR,3,'omitnan');
semLDR = std(LDR,[],3)./sqrt(size(LDR,3));
for u = 1:size(LDR,1)
    pvalLDR(u,1) = ranksum(squeeze(LDR(u,1,:)),squeeze(LDR(u,2,:)));

    mV = meanvec(squeeze(ANG_DIR(u,1,:)),ones(size(ANG_DIR,3),1));
    mANG_DIR(u,1) = mV.angDir;
    dDir = abs(squeeze(ANG_DIR(u,1,:)) - mANG_DIR(u,1)); dDir(dDir>180) = 360-dDir(dDir>180);
    semANG_DIR(u,1) = sqrt(sum(dDir.^2)/(nReps-1))/sqrt(nReps);

    mV = meanvec(squeeze(ANG_DIR(u,2,:)),ones(size(ANG_DIR,3),1));
    mANG_DIR(u,2) = mV.angDir;
    dDir = abs(squeeze(ANG_DIR(u,2,:)) - mANG_DIR(u,2)); dDir(dDir>180) = 360-dDir(dDir>180);
    semANG_DIR(u,2) = sqrt(sum(dDir.^2)/(nReps-1))/sqrt(nReps);
end
% mANG_DIR = mean(ANG_DIR,3,'omitnan');
% semANG_DIR = std(ANG_DIR,[],3)./sqrt(size(ANG_DIR,3));

mLDR_s = mean(LDR_s,3,'omitnan');
semLDR_s = std(LDR_s,[],3)./sqrt(size(LDR_s,3));
for u = 1:size(LDR_s,1)
    pvalLDR_s(u,1) = ranksum(squeeze(LDR_s(u,1,:)),squeeze(LDR_s(u,2,:)));
end
% mANG_DIR_s = mean(ANG_DIR_s,3,'omitnan');
% semANG_DIR_s = std(ANG_DIR_s,[],3)./sqrt(size(ANG_DIR_s,3));

mLOR = mean(LOR,3,'omitnan');
semLOR = std(LOR,[],3)./sqrt(size(LOR,3));
for u = 1:size(LOR,1)
    pvalLOR(u,1) = ranksum(squeeze(LOR(u,1,:)),squeeze(LOR(u,2,:)));
end
% mANG_ORI = mean(ANG_ORI,3,'omitnan');
% semANG_ORI = std(ANG_ORI,[],3)./sqrt(size(ANG_ORI,3));

mLOR_s = mean(LOR_s,3,'omitnan');
semLOR_s = std(LOR_s,[],3)./sqrt(size(LOR_s,3));
for u = 1:size(LOR_s,1)
    pvalLOR_s(u,1) = ranksum(squeeze(LOR_s(u,1,:)),squeeze(LOR_s(u,2,:)));
end
% mANG_ORI_s = mean();

mDSI = mean(DSI,3,'omitnan');
semDSI = std(DSI,[],3)./sqrt(size(DSI,3));
for u = 1:size(DSI,1)
    pvalDSI(u,1) = ranksum(squeeze(DSI(u,1,:)),squeeze(DSI(u,2,:)));
end

mDSI_sRep = mean(DSI_sRep,3,'omitnan');
semDSI_sRep = std(DSI_sRep,[],3)./sqrt(size(DSI_sRep,3));
for u = 1:size(DSI_sRep,1)
    pvalDSI_sRep(u,1) = ranksum(squeeze(DSI_sRep(u,1,:)),squeeze(DSI_sRep(u,2,:)));
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

mSize = 7;
agShapes = {'+','square','diamond','^'};

exUs = [4 25];
exUclrs = {[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880]};
for u = 10

    %Tuning Curve
    x = conds;
    y1 = tc(u,:,1);
    sem1 = tc_sem(u,:,1);
    y2 = tc(u,:,2);
    sem2 = tc_sem(u,:,2);

%     figure; hold on
%     errorbar(x,y1,sem1,'k','LineWidth',1)
%     errorbar(x,y2,sem2,'c','LineWidth',1)
%     plot(conds,y2*(rPref(u,1)/rPref(u,2)),'c--','LineWidth',1)
%     % legend('cntrl','cool','scaled cool')
%     xticks([0 90 180 270])
%     xlim([0 360])
%     xlabel('dir of motion (deg)')
%     ylabel('firing rate')

    %polar TC
    figure;hold on
    %polar axes
    rtic = 0;
    while rtic<max([y1,y2])
        rtic = rtic+5;
        [axT,axR] = pol2cart(deg2rad(0:360),ones(1,361)*rtic);
        plot(axT,axR,'k-')
        rticMax = rtic;
    end
%     rtic = 0;
%     while rtic<rticMax
%         rtic = rtic+1;
%         [axT,axR] = pol2cart(deg2rad(0:360),ones(1,361)*rtic);
%         plot(axT,axR,'k:')
%     end
%     ttic = 0;
%     while ttic<360
%         [axT,axR] = pol2cart(deg2rad([ttic ttic]),[0 rticMax]);
%         plot(axT,axR,'k:')
%         ttic = ttic+22.5;
%     end
    ttic = 0;
    while ttic<360
        [axT,axR] = pol2cart(deg2rad([ttic ttic]),[0 rticMax]);
        plot(axT,axR,'k--')
        ttic = ttic+45;
    end

    mv1 = meanvec(x,y1);
    mv2 = meanvec(x,y2);
    [mv1_x,mv1_y] = pol2cart(repmat(deg2rad(mv1.angDir),1,2),[0 rticMax*mv1.ldr]);
    [mv2_x,mv2_y] = pol2cart(repmat(deg2rad(mv2.angDir),1,2),[0 rticMax*mv2.ldr]);
    
    y1 = [y1 y1(1)]; sem1 = [sem1 sem1(1)]; y2 = [y2 y2(1)];sem2 = [sem2 sem2(1)] ; x = [x x(1)];
    tmp1 = [y1+sem1 fliplr(y1-sem1)];tmp1(tmp1<0) = 0;
    [semX1,semY1] = pol2cart(deg2rad([x fliplr(x)]),tmp1);
    tmp2 = [y2+sem2 fliplr(y2-sem2)];tmp2(tmp2<0) = 0;
    [semX2,semY2] = pol2cart(deg2rad([x fliplr(x)]),tmp2);
    patch(semX1,semY1,'k','EdgeColor','none','FaceAlpha',0.2)
    patch(semX2,semY2,'c','EdgeColor','none','FaceAlpha',0.2)
    plot(semX1,semY1,'k')
    plot(semX2,semY2,'c')
    [x2_scl,y2_scl] = pol2cart(deg2rad([x x(1)]),[y2 y2(1)]*(rPref(u,1)/rPref(u,2)));
    [x1,y1] = pol2cart(deg2rad([x x(1)]),[y1 y1(1)]);
    [x2,y2] = pol2cart(deg2rad([x x(1)]),[y2 y2(1)]);
    plot(x1,y1,'k','LineWidth',2)
    plot(x2,y2,'c','LineWidth',2)
    plot(x2_scl,y2_scl,'c--','LineWidth',2)
    plot(mv1_x,mv1_y,'k','LineWidth',3);
    plot(mv2_x,mv2_y,'c','LineWidth',3);

    box on
    axis square
    title([dat{1}.exptName{u}(1:5) ' ' dat{1}.area{u} ' ' dat{1}.uInfo{u} num2str(dat{1}.uID(u))])
%     saveas(gcf,fullfile(dataFold,'Figures',proj,'matchedSU',['su' num2str(u)]),'fig')
%     saveas(gcf,fullfile(dataFold,'Figures',proj,'matchedSU',['su' num2str(u)]),'svg')

    
%     figure; hold on
%     countTrials = 0;
%     for i = 1:3
%         if i==1
%             clr = 'k';
%         elseif i==2
%             clr = 'c';
%         elseif i==3
%             clr = 'k';
%         end
%         x = dat{i}.spkTimes{u}(1,:);
%         y = dat{i}.spkTimes{u}(2,:);
%         y = y+countTrials;
%         plot(x,y,[clr '.'],'MarkerSize',4)
%         countTrials = countTrials+max(dat{i}.fr(u).trialNum,[],'all');
%     end
%     xlim([-1 2])
%     xlabel('time (sec)')
%     ylim([0 countTrials+1])
%     ylabel('trial number')
%     patch([0 1 1 0],[0 0 countTrials+1 countTrials+1],'k','EdgeColor','none','FaceAlpha',0.2)
%     box on
%     axis square

end

figure; hold on
% subplot(3,3,4); hold on
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};

% x = [ldr(:,1),mLDR(:,1)]';
% y = [ldr(:,2),mLDR(:,2)]';
% dLDR(1,:) = sqrt(sum([diff(x).^2;diff(y).^2]));
% % plot(x,y,'-','Color',colors{1},'LineWidth',2)
% 
% x = [ldr_s(:,1),mLDR_s(:,1)]';
% y = [ldr_s(:,2),mLDR_s(:,2)]';
% dLDR(2,:) = sqrt(sum([diff(x).^2;diff(y).^2]));
% % plot(x,y,'--','Color',colors{2},'LineWidth',2)
% 
% x = [mLDR_s(:,1),mLDR(:,1)]';
% y = [mLDR_s(:,2),mLDR(:,2)]';
% dLDR(3,:) = sqrt(sum([diff(x).^2;diff(y).^2]));
% % plot(x,y,'-.','Color',colors{3},'LineWidth',2)
% 
% x = [ldr(:,1),ldr_s(:,1)]';
% y = [ldr(:,2),ldr_s(:,2)]';
% dLDR(4,:) = sqrt(sum([diff(x).^2;diff(y).^2]));
% % plot(x,y,':','Color',colors{4},'LineWidth',2)
% 
% % cdfP = cdfplot(dLDR(1,:));
% % cdfP.LineStyle = '-';
% % % cdfP.Color = 'k';
% % cdfP.LineWidth = 3;
% % cdfP = cdfplot(dLDR(2,:));
% % cdfP.LineStyle = '--';
% % % cdfP.Color = 'k';
% % cdfP.LineWidth = 3;
% % cdfP = cdfplot(dLDR(3,:));
% % cdfP.LineStyle = '-.';
% % % cdfP.Color = 'k';
% % cdfP.LineWidth = 3;
% % cdfP = cdfplot(dLDR(4,:));
% % cdfP.LineStyle = ':';
% % % cdfP.Color = 'k';
% % cdfP.LineWidth = 3;
% % axis square
% % box on
% % xlabel('dL_D_i_r (euc. dist.)')
% % xlim([0 0.2])
% % ylabel('percentile')

% clr = 'b';
% x = ldr(:,1);
% y = ldr(:,2);
% plot(x,y,[clr 'o'],'MarkerSize',mSize,'MarkerFaceColor',clr);
% 
% clr = 'c';
% x = mLDR(:,1);
% errX = semLDR(:,1);
% y = mLDR(:,2);
% errY = semLDR(:,2);
% sigIdx = pvalLDR<0.05;
% p(1) = errorbar(x(sigIdx),y(sigIdx),errY(sigIdx),errY(sigIdx),errX(sigIdx),errX(sigIdx),[clr 'o'],'MarkerFaceColor',clr,'CapSize',0);
% p(1).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unit',find(sigIdx));
% p(2) = errorbar(x(~sigIdx),y(~sigIdx),errY(~sigIdx),errY(~sigIdx),errX(~sigIdx),errX(~sigIdx),[clr 'o'],'MarkerFaceColor','w','CapSize',0);
% p(2).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unit',find(~sigIdx));
% 
% clr = 'r';
% x = ldr_s(:,1);
% y = ldr_s(:,2);
% plot(x,y,[clr 'o'],'MarkerSize',mSize,'MarkerFaceColor',clr);

clr = 'k';
x = mLDR_s(:,1);
errX = semLDR_s(:,1);
y = mLDR_s(:,2);
errY = semLDR_s(:,2);
sigIdx = pvalLDR_s<0.05;
p(1) = errorbar(x(sigIdx),y(sigIdx),errY(sigIdx),errY(sigIdx),errX(sigIdx),errX(sigIdx),[clr 'o'],'MarkerFaceColor',clr,'CapSize',0);
p(1).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unit',find(sigIdx));
p(2) = errorbar(x(~sigIdx),y(~sigIdx),errY(~sigIdx),errY(~sigIdx),errX(~sigIdx),errX(~sigIdx),[clr 'o'],'MarkerFaceColor','w','CapSize',0);
p(2).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unit',find(~sigIdx));

% plot(x(exUs(1)),y(exUs(1)),'o','MarkerSize',mSize,'Color',exUclrs{1},'LineWidth',1)
% plot(x(exUs(2)),y(exUs(2)),'o','MarkerSize',mSize,'Color',exUclrs{2},'LineWidth',1)

plot([0 1],[0 1],'k--')
xlabel('mean vector length, pre-cooling')
ylabel('mean vector length, cooling')
title('L_D_i_r')
% legend(p,{'p < 0.05','p > 0.05'},'Location','southeast')
box on
axis square



figure; hold on
% % subplot(3,3,5); hold on

x = [dsi(:,1),mDSI(:,1)]';
y = [dsi(:,2),mDSI(:,2)]';
dDSI(1,:) = sqrt(sum([diff(x).^2;diff(y).^2]));
% plot(x,y,'-','Color',colors{1},'LineWidth',2)

x = [dsi_s(:,1),mDSI_sRep(:,1)]';
y = [dsi_s(:,2),mDSI_sRep(:,2)]';
dDSI(2,:) = sqrt(sum([diff(x).^2;diff(y).^2]));
% plot(x,y,'--','Color',colors{2},'LineWidth',2)

x = [mDSI_sRep(:,1),mDSI(:,1)]';
y = [mDSI_sRep(:,2),mDSI(:,2)]';
dDSI(3,:) = sqrt(sum([diff(x).^2;diff(y).^2]));
% plot(x,y,'-.','Color',colors{3},'LineWidth',2)

x = [dsi(:,1),dsi_s(:,1)]';
y = [dsi(:,2),dsi_s(:,2)]';
dDSI(4,:) = sqrt(sum([diff(x).^2;diff(y).^2]));
% plot(x,y,':','Color',colors{4},'LineWidth',2)

% cdfP = cdfplot(dDSI(1,:));
% cdfP.LineStyle = '-';
% % cdfP.Color = 'k';
% cdfP.LineWidth = 3;
% cdfP = cdfplot(dDSI(2,:));
% cdfP.LineStyle = '--';
% % cdfP.Color = 'k';
% cdfP.LineWidth = 3;
% cdfP = cdfplot(dDSI(3,:));
% cdfP.LineStyle = '-.';
% % cdfP.Color = 'k';
% cdfP.LineWidth = 3;
% cdfP = cdfplot(dDSI(4,:));
% cdfP.LineStyle = ':';
% % cdfP.Color = 'k';
% cdfP.LineWidth = 3;
% axis square
% box on
% xlabel('dDSI (euc. dist.)')
% ylabel('percentile')

% clr = 'b';
% x = dsi(:,1);
% y = dsi(:,2);
% plot(x,y,[clr 'o'],'MarkerSize',mSize,'MarkerFaceColor',clr)
% 
% clr = 'c';
% x = mDSI(:,1);
% errX = semDSI(:,1);
% y = mDSI(:,2);
% errY = semDSI(:,2);
% sigIdx = pvalDSI<0.05;
% p(1) = errorbar(x(sigIdx),y(sigIdx),errY(sigIdx),errY(sigIdx),errX(sigIdx),errX(sigIdx),[clr 'o'],'MarkerFaceColor',clr,'CapSize',0);
% p(1).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unit',find(sigIdx));
% p(2) = errorbar(x(~sigIdx),y(~sigIdx),errY(~sigIdx),errY(~sigIdx),errX(~sigIdx),errX(~sigIdx),[clr 'o'],'MarkerFaceColor','w','CapSize',0);
% p(2).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unit',find(~sigIdx));
% 
% clr = 'r';
% x = dsi_s(:,1);
% y = dsi_s(:,2);
% plot(x,y,[clr 'o'],'MarkerSize',mSize,'MarkerFaceColor',clr)

clr = 'k';
x = mDSI_sRep(:,1);
errX = semDSI_sRep(:,1);
y = mDSI_sRep(:,2);
errY = semDSI_sRep(:,2);
sigIdx = pvalDSI_sRep<0.05;
p(1) = errorbar(x(sigIdx),y(sigIdx),errY(sigIdx),errY(sigIdx),errX(sigIdx),errX(sigIdx),[clr 'o'],'MarkerFaceColor',clr,'CapSize',0);
p(1).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unit',find(sigIdx));
p(2) = errorbar(x(~sigIdx),y(~sigIdx),errY(~sigIdx),errY(~sigIdx),errX(~sigIdx),errX(~sigIdx),[clr 'o'],'MarkerFaceColor','w','CapSize',0);
p(2).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unit',find(~sigIdx));

% plot(x(exUs(1)),y(exUs(1)),'o','MarkerSize',mSize,'Color',exUclrs{1},'LineWidth',1)
% plot(x(exUs(2)),y(exUs(2)),'o','MarkerSize',mSize,'Color',exUclrs{2},'LineWidth',1)

plot([0 1],[0 1],'k--')
xlim([0 1]);ylim([0 1])
xlabel('DSI, pre-cooling')
ylabel('DSI, cooling')
title('DSI')
% legend(p,{'p < 0.05','p > 0.05'},'Location','southeast')
box on
axis square



% subplot(3,3,6); hold on
figure; hold on
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


% subplot(3,2,5); hold on
figure; hold on
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

% subplot(3,2,6); hold on
figure; hold on
boxchart(categorical(d.ori),d.si,'notch','on','BoxFaceColor','k')
ylabel('SI')
% boxchart(categorical(d.ori),d.siL,'notch','on')
% ylabel('SI (log)')
yline(0,'k--')
xlabel('Angle (deg. relative to preferred)')
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
