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
%     if ~strcmp(animals{a},'febh5')
%         continue
%     end
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
%         splitIntan(fullfile(dataFold,'Ephys'),animals{a},mergeId,probeId,'BRN')
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
        
        pvalANOVA(u,e) = anova1(rMat,[],"off");

        %tuning curve
        rMean = mean(rMat,1,'omitnan');
        rSem = std(rMat,[],1,'omitnan')/sqrt(size(rMat,1));
        pvalTUNE(u,e) = anova1(rMat,[],'off');

        wHan=hann(3);
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
            
            [~,BW(u,e,rep)] = compute_bw(conds,rRep);

            % [BW(u,e,rep),BW_s(u,e,rep)] = compute_bw(conds,rRep);


            wHan=hann(3);
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
            [c_aligned,tc_aligned(u,:,e),alignIdx(u,:,e)] = alignDirTuning(conds,rMean);
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
semLDR = std(LDR,[],3)/sqrt(nReps);
for u = 1:size(LDR,1)
    pvalLDR(u,1) = ranksum(squeeze(LDR(u,1,:)),squeeze(LDR(u,2,:)));

    mV = meanvec(squeeze(ANG_DIR(u,1,:)),ones(size(ANG_DIR,3),1));
    mV_s = meanvec(squeeze(ANG_DIR_s(u,1,:)),ones(size(ANG_DIR_s,3),1));
    mANG_DIR(u,1) = mV.angDir;
    mANG_DIR_s(u,1) = mV_s.angDir;
    dDir = abs(squeeze(ANG_DIR(u,1,:)) - mANG_DIR(u,1)); dDir(dDir>180) = 360-dDir(dDir>180);
    semANG_DIR(u,1) = sqrt(sum(dDir.^2)/(nReps-1))/sqrt(nReps);

    mV = meanvec(squeeze(ANG_DIR(u,2,:)),ones(size(ANG_DIR,3),1));
    mV_s = meanvec(squeeze(ANG_DIR_s(u,2,:)),ones(size(ANG_DIR_s,3),1));
    mANG_DIR(u,2) = mV.angDir;
    mANG_DIR_s(u,2) = mV_s.angDir;
    dDir = abs(squeeze(ANG_DIR(u,2,:)) - mANG_DIR(u,2)); dDir(dDir>180) = 360-dDir(dDir>180);
    semANG_DIR(u,2) = sqrt(sum(dDir.^2)/(nReps-1))/sqrt(nReps);
end

mLDR_s = mean(LDR_s,3,'omitnan');
semLDR_s = std(LDR_s,[],3)/sqrt(nReps);
for u = 1:size(LDR_s,1)
    pvalLDR_s(u,1) = ranksum(squeeze(LDR_s(u,1,:)),squeeze(LDR_s(u,2,:)));
end

mLOR = mean(LOR,3,'omitnan');
semLOR = std(LOR,[],3)./sqrt(nReps);
for u = 1:size(LOR,1)
    pvalLOR(u,1) = ranksum(squeeze(LOR(u,1,:)),squeeze(LOR(u,2,:)));
end

mLOR_s = mean(LOR_s,3,'omitnan');
semLOR_s = std(LOR_s,[],3)/sqrt(nReps);
for u = 1:size(LOR_s,1)
    pvalLOR_s(u,1) = ranksum(squeeze(LOR_s(u,1,:)),squeeze(LOR_s(u,2,:)));
end

mBW = mean(BW,3,'omitnan');
semBW = std(BW,[],3)/sqrt(nReps);
for u = 1:size(BW,1)
    pvalBW(u,1) = ranksum(squeeze(BW(u,1,:)),squeeze(BW(u,2,:)));
end

mDSI = mean(DSI,3,'omitnan');
semDSI = std(DSI,[],3)./sqrt(nReps);
for u = 1:size(DSI,1)
    pvalDSI(u,1) = ranksum(squeeze(DSI(u,1,:)),squeeze(DSI(u,2,:)));
end

mDSI_s = mean(DSI_sRep,3,'omitnan');
semDSI_sRep = std(DSI_sRep,[],3)./sqrt(nReps);
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
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};

% plot example neurons
plr=1;
exUs = [4 10 5 23];
% exUs = [10 47];
% exUs = 1:nU;
for u = exUs

    figure;hold on
    for e = 1:2
        if e == 2 
            clr = 'c';
        else
            clr = 'k';
        end
        %Tuning Curve
        x = conds;
        y = tc(u,:,e);
        sclF = (rPref(u,1)/rPref(u,2));
        mv_x = mANG_DIR(u,e);
        mv_y = mLDR(u,e);
        err = tc_sem(u,:,e);
        tmpY = [y+err y(1)+err(1) y(1)-err(1) fliplr(y-err)];
        tmpY(tmpY<0) = 0;
        if plr == 1
            % construct polar axes
            rtic = 0;
            while rtic<max([tc(u,:,1),tc(u,:,2)])
                rtic = rtic+5;
                [axT,axR] = pol2cart(deg2rad(0:360),ones(1,361)*rtic);
                plot(axT,axR,'k-')
                rticMax = rtic;
            end
            ttic = 0;
            while ttic<360
                [axT,axR] = pol2cart(deg2rad([ttic ttic]),[0 rticMax]);
                plot(axT,axR,'k--')
                ttic = ttic+45;
            end
            % convert to polar coordinates
            [mv_xP,mv_yP] = pol2cart(repmat(deg2rad(mv_x),1,2),[0 rticMax*mv_y]);
            [err_xP,err_yP] = pol2cart(deg2rad([x x(1) x(1) fliplr(x)]),tmpY);
            [xP,yP] = pol2cart(deg2rad(x),y);
            % connect ends
            xP = [xP xP(1)];
            yP = [yP yP(1)];
            % plot
            patch(err_xP,err_yP,clr,'EdgeColor','none','FaceAlpha',0.2)
            plot(err_xP,err_yP,clr)
            plot(xP,yP,clr,'LineWidth',2)
            plot(mv_xP,mv_yP,clr,'LineWidth',3)
            % if cooling condition, plot scaled tuning curve
            if e == 2
                [scldTC_xP,scldTC_yP] = pol2cart(deg2rad(x),y*sclF);
                scldTC_xP = [scldTC_xP scldTC_xP(1)];
                scldTC_yP = [scldTC_yP scldTC_yP(1)];
                plot(scldTC_xP,scldTC_yP,'c--','LineWidth',2)
            end
        else
            errorbar([x x(1)+360],[y y(1)],[err err(1)],clr,'LineWidth',1)
            if e == 2
                plot([x x(1)+360],[y y(1)]*sclF,[clr '--'],'LineWidth',1)
            end
            % legend('cntrl','cool','scaled cool')
            xticks([0 90 180 270])
            xlim([0 360])
            xlabel('dir of motion (deg)')
            ylabel('firing rate')
        end
        box on
        axis square
    end

    box on
    axis square
    title([dat{1}.exptName{u}(1:5) ' ' dat{1}.area{u} ' ' dat{1}.uInfo{u} num2str(dat{1}.uID(u))])
%     saveas(gcf,fullfile(dataFold,'Figures',proj,'matchedSU',['su' num2str(u)]),'fig')
%     saveas(gcf,fullfile(dataFold,'Figures',proj,'matchedSU',['su' num2str(u)]),'svg')

end

% scatter plot - Ldir
figure; hold on
clr = 'k';
x = mLDR(:,1);
errX = semLDR(:,1);
y = mLDR(:,2);
errY = semLDR(:,2);
sigIdx = pvalLDR<0.05;
p(2) = errorbar(x(~sigIdx),y(~sigIdx),errY(~sigIdx),errY(~sigIdx),errX(~sigIdx),errX(~sigIdx),[clr 'o'],'MarkerFaceColor','w','CapSize',0);
p(2).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unit',find(~sigIdx));
p(1) = errorbar(x(sigIdx),y(sigIdx),errY(sigIdx),errY(sigIdx),errX(sigIdx),errX(sigIdx),[clr 'o'],'MarkerFaceColor',clr,'CapSize',0);
p(1).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unit',find(sigIdx));
plot([0 1],[0 1],'k--')
xlabel('mean vector length, pre-cooling')
ylabel('mean vector length, cooling')
title('L_D_i_r')
legend(p,{'p < 0.05','p > 0.05'},'Location','southeast')
box on
axis square

% scatter plot - DSI
figure; hold on
clr = 'k';
x = mDSI_s(:,1);
errX = semDSI_sRep(:,1);
y = mDSI_s(:,2);
errY = semDSI_sRep(:,2);
sigIdx = pvalDSI_sRep<0.05;
p(2) = errorbar(x(~sigIdx),y(~sigIdx),errY(~sigIdx),errY(~sigIdx),errX(~sigIdx),errX(~sigIdx),[clr 'o'],'MarkerFaceColor','w','CapSize',0);
p(2).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unit',find(~sigIdx));
p(1) = errorbar(x(sigIdx),y(sigIdx),errY(sigIdx),errY(sigIdx),errX(sigIdx),errX(sigIdx),[clr 'o'],'MarkerFaceColor',clr,'CapSize',0);
p(1).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unit',find(sigIdx));
plot([0 1],[0 1],'k--')
xlim([0 1]);ylim([0 1])
xlabel('DSI, pre-cooling')
ylabel('DSI, cooling')
title('DSI')
% legend(p,{'p < 0.05','p > 0.05'},'Location','southeast')
box on
axis square

% polar histogram - change in pref dir
figure;
dDir = mANG_DIR(:,1)-mANG_DIR(:,2);
dDir(dDir>180) = dDir(dDir>180)-360;
dDir(dDir<-180) = dDir(dDir<-180)+360;
sigIdx = pvalANOVA(:,1)<0.05 & pvalANOVA(:,2)<0.05;
% histogram(dDir,-180:10:180)
% xticks([-180:90:180])
% box on
% axis square
binW = 22.5;
bins = -binW/2:binW:360-(binW/2);
polarhistogram(deg2rad(dDir),deg2rad(bins));hold on
p(2) = polarscatter(deg2rad(dDir(~sigIdx)),ones(size(dDir(~sigIdx)))*25,'ko','MarkerFaceColor','w');
p(2).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unit',find(~sigIdx));
p(1) = polarscatter(deg2rad(dDir(sigIdx)),ones(size(dDir(sigIdx)))*25,'ko','MarkerFaceColor','k');
p(1).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unit',find(sigIdx));
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.ThetaTick = bins;
ax.RTick = [0:5:25];
% ax.ThetaTickLabel = {'0','45','90','135','+/-180','-135','-90','-45'};
ax.GridAlpha = 1;
title('Difference in preferred direction')

% scatter plot - bandwidth
clr = 'k';
figure; hold on
plot([0 100],[0 100],'k--')
x = mBW(:,1);
errX = semBW(:,1);
y = mBW(:,2);
errY = semBW(:,2);
sigIdx = pvalBW<0.05;
p(2) = errorbar(x(~sigIdx),y(~sigIdx),errY(~sigIdx),errY(~sigIdx),errX(~sigIdx),errX(~sigIdx),[clr 'o'],'MarkerFaceColor','w','CapSize',0);
p(2).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unit',find(~sigIdx));
p(1) = errorbar(x(sigIdx),y(sigIdx),errY(sigIdx),errY(sigIdx),errX(sigIdx),errX(sigIdx),[clr 'o'],'MarkerFaceColor',clr,'CapSize',0);
p(1).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unit',find(sigIdx));
xlabel('BW, pre-cooling')
ylabel('BW, cooling')
title('Bandwidth')
legend(p,{'p < 0.05','p > 0.05'},'Location','southeast')
box on
axis square

% box plot - Rpref per direction/manipulation
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

% box plot - SI per direction
figure; hold on
boxchart(categorical(d.ori),d.si,'notch','on','BoxFaceColor','k')
ylabel('SI')
% boxchart(categorical(d.ori),d.siL,'notch','on')
% ylabel('SI (log)')
for u = exUs 
    plot((r2(u,:)-r1(u,:))./(r2(u,:)+r1(u,:)),'-o');
end
yline(0,'k--')
xlabel('Angle (deg. relative to preferred)')
box on



r1 = tc(:,:,1);
r2 = tc(:,:,2);
si = (r2-r1)./(r2+r1);
for u = 1:size(si,1)
    if ismember(u,find(pvalANOVA(:,1)<0.05))
        [cAl,~,idx] = alignDirTuning(conds,r1(u,:));
        siAl(u,:) = si(u,idx);
        r1Al(u,:) = r1(u,idx);
        r2Al(u,:) = r2(u,idx);
%     elseif ismember(u,find(pvalANOVA(:,2)<0.05))
%         [cAl,~,idx] = alignDirTuning(conds,r2(u,:));
%         siAl(u,:) = si(u,idx);
%         r1Al(u,:) = r1(u,idx);
%         r2Al(u,:) = r2(u,idx);
    else
        siAl(u,:) = nan(1,size(si,2)+1);
        r1Al(u,:) = nan(1,size(si,2)+1);
        r2Al(u,:) = nan(1,size(si,2)+1);
    end

    ad = unique(abs(cAl));
    for c = 1:length(ad)
        siAlH(u,c) = mean(siAl(u,abs(cAl)==ad(c)),'omitnan');
    end


end

figure; hold on
boxchart(siAlH)

% figure;hold on
% colororder({'r'})
% u = 8;
% plot(r1AL(u,:),'k','LineWidth',2)
% plot(r2AL(u,:),'c','LineWidth',2)
% yyaxis right
% plot(siAL(u,:),'r','LineWidth',2)
% yline(0,'r--')

%% OLD

% % scatter plot - bandwidth (old)
% figure; hold on
% x = bwS(:,2);
% y = bwS(:,1);
% lims = [10 60];
% out = x>lims(2)|y>lims(2);
% % plot(x,y,'k.','MarkerSize',4)
% for ag = 1:length(ageGroups)
%     plot(x(uAG==ag),y(uAG==ag),['k' agShapes{ag}],'MarkerSize',mSize,'LineWidth',1)
% end
% x(x>lims(2)) = lims(2);
% y(y>lims(2)) = lims(2);
% for ag = 1:length(ageGroups)
%     plot(x(uAG==ag & out'),y(uAG==ag & out'),['r' agShapes{ag}],'MarkerSize',mSize,'LineWidth',1)
% end
% plot(x(exUs(1)),y(exUs(1)),'o','MarkerSize',mSize,'Color',exUclrs{1},'LineWidth',1)
% plot(x(exUs(2)),y(exUs(2)),'o','MarkerSize',mSize,'Color',exUclrs{2},'LineWidth',1)
% plot([0 200],[0 200],'k--')
% xlabel('V1 cooled') 
% xlim(lims)
% ylim(lims)
% ylabel('control')
% title('bandwidth (smooth)')
% box on
% axis square