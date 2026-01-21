clear all
close all

%% Load data

anaMode = 'SU';

% dataFold = '/Volumes/Lab drive/Brandon/data/dataSets/DSdev';
% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/DSdev';
dataFold = '/Volumes/NielsenHome2/Brandon/data/dataSets/DSdev';
% dataFold = 'F:\Brandon\data\dataSets\DSdev';
% dataFold = 'Y:\Brandon\data\dataSets\DSdev';
load(fullfile(dataFold,['DSdev_' anaMode 'dataSet.mat']))
dir = load(fullfile(dataFold,'anaRSA_dir.mat'));
ori = load(fullfile(dataFold,'anaRSA_ori.mat'));

%% Organize data

areas = {'V1','PSS'};
ageGroups = {[29 32],[33 36],[37 max(projTbl.age)]};
% ageGroups = {[28 32],[29 33],[30 34],[31 35],[32 36],[33 37],[34 38],[35 39],[36 40],[37 41],[38 42],[39 43],[40 44],[41 300]};

nAG = length(ageGroups);
nAR = length(areas);
Rtime = cell(nAR,nAG);
Rtrial = cell(nAR,nAG);
binSpacing = 0.05;
bins = -1:binSpacing:2;
binOverlap = 0.025;
binSize = binSpacing+(2*binOverlap);
% binRight = bins+(binSize/2);
% binLeft = bins-(binSize/2);
binRight = bins+binSize;
binLeft = bins;
nBins = length(bins);
for ar = 1:nAR
for ag = 1:nAG

    ageLims = ageGroups{ag};
    areaIdx = strcmp(projTbl.recSite,areas{ar});
    ageLimIdx = projTbl.age>=ageLims(1) & projTbl.age<=ageLims(2);
    dat{ar,ag} = vertcat(projTbl(areaIdx & ageLimIdx,:).sumStats{:});

    %only take units that pass inclusion criteria (screenUnits)
    dat{ar,ag} = dat{ar,ag}(screenUnits(dat{ar,ag},anaMode),:);
    
    %only take units with ori12
    for u = 1:height(dat{ar,ag})
        isOri12(u) = size(dat{ar,ag}.condition{u},2)==12;
        isOri16(u) = size(dat{ar,ag}.condition{u},2)==16;
    end
    dat{ar,ag} = dat{ar,ag}(isOri12,:);
    clear isOri12 isOri16

    nU(ar,ag) = height(dat{ar,ag});

    %sort units by direction preference
    [~,sortIdx] = sort(dat{ar,ag}.oriPref);
    dat{ar,ag} = dat{ar,ag}(sortIdx,:);

    %make response matrix
    for u = 1:height(dat{ar,ag})

        hasBlank = ~isempty(dat{ar,ag}.rBlank{u});
        oriIdx = contains(dat{ar,ag}.paramKey{u},'ori');
        sizeIdx = contains(dat{ar,ag}.paramKey{u},'size');
        posIdx = contains(dat{ar,ag}.paramKey{u},'pos');
        hasSize = sum(sizeIdx|posIdx)>0;
        if hasSize
            test=1;
            cndInclude = find(dat{ar,ag}.cndKey{u}(:,(sizeIdx|posIdx))>100);
        else
            cndInclude = 1:size(dat{ar,ag}.cndKey{u},1);
        end
%         cndInclude = cndInclude(ismember(dat{ar,ag}.cndKey{u}(cndInclude,oriIdx),[0 90 180 270]));

        conds = dat{ar,ag}.cndKey{u}(cndInclude,oriIdx);
        nConds = length(conds);
        nReps = 5;
        nTrials = nReps*nConds;

        spks = dat{ar,ag}.spkTimes{u}(1,:);
        spkTrial = dat{ar,ag}.spkTimes{u}(2,:);
        trialIds = dat{ar,ag}.fr(u).trialNum(1:nReps,cndInclude);
        trialIds = trialIds(:);
        trialConds = dat{ar,ag}.fr(u).trialCond(1:nReps,cndInclude);
        trialConds = trialConds(:);
        if min(trialConds)~=1
            trialConds = trialConds-(min(trialConds)-1);
        end

        dirs = conds;
        Cdir = reshape(repmat(dirs,1,nReps)',nReps*nConds,1);%dirs(trialConds);
        oris = unique(mod(dirs,180));
        Cori = mod(Cdir,180);

        %trial wise stim response
        rTmp = dat{ar,ag}.response{u}(1:nReps,(cndInclude-(min(cndInclude)-1)));
        rTmp(rTmp<0) = 0;
        Rtrial{ar,ag}(:,u) = rTmp(:);

        %trial averaged response (by direction class)
        for d = 1:length(unique(dirs))
            RmeanDir{ar,ag}(d,u) = mean( Rtrial{ar,ag}(Cdir==dirs(d),u) );
        end

        %trial averaged response (by orientation class)
        for o = 1:length(unique(oris))
            RmeanOri{ar,ag}(o,u) = mean( Rtrial{ar,ag}(Cori==oris(o),u) );
        end

        %trial wise time resolved
        for t = 1:nTrials
            tId = trialIds(t);
%             psth = histcounts(spks(spkTrial==tId),bins);
%             psth = psth/binSize;
            for b = 1:nBins
                psth(b) = sum(spks(spkTrial==tId)>binLeft(b)&spks(spkTrial==tId)<=binRight(b))/binSize;
            end
            Rtime{ar,ag}(t,u,:) = psth;
            clear psth
        end

        for d = 1:length(unique(dirs))
            Rtime_meanDir{ar,ag}(d,u,:) = mean( squeeze(Rtime{ar,ag}(Cdir==dirs(d),u,:)) );
        end

        for o = 1:length(unique(oris))
            Rtime_meanOri{ar,ag}(o,u,:) = mean( squeeze(Rtime{ar,ag}(Cori==oris(o),u,:)) );
        end



    end
        
end
end

ar = 1; ag = 3;


figure;hold on
x = Rtrial{ar,ag};
y = Cdir;
uY = unique(y);
imagesc(x)
xlabel('neuron')
ylabel('trial')
axis square
axis tight
title([areas{ar} ' P' num2str(ageGroups{ag}(1)) ' - ' num2str(ageGroups{ag}(2))])

figure;hold on
[coeff,score] = pca(x);
clrs = hsv(length(uY));
for c = 1:length(uY)
    idx = y==uY(c);
    idx2 = dirs == dirs(c);
    plot3(score(idx,1),score(idx,2),score(idx,3),'o','MarkerFaceColor',clrs(c,:))
end
plot3(score(:,1),score(:,2),score(:,3),'ko')
xlabel('PC1');ylabel('PC2');zlabel('PC3')
axis square
box on
clear score
title([areas{ar} ' P' num2str(ageGroups{ag}(1)) ' - ' num2str(ageGroups{ag}(2))])


figure;hold on
x = zscore(RmeanDir{ar,ag});
y = dirs;
uY = unique(y);
[coeff,score] = pca(x);
clrs = hsv(length(uY));
for c = 1:length(uY)
    idx = y==uY(c);
    plot3(score(idx,1),score(idx,2),score(idx,3),'o','MarkerFaceColor',clrs(c,:))
end
plot3([score(:,1);score(1,1)],[score(:,2);score(1,2)],[score(:,3);score(1,3)],'k-o')
xlabel('PC1');ylabel('PC2');zlabel('PC3')
axis square
box on
clear score
title([areas{ar} ' P' num2str(ageGroups{ag}(1)) ' - ' num2str(ageGroups{ag}(2))])


figure('Position',[100 100 1000 600])
for ar = 1:nAR
    for ag = 1:nAG
        subplot(nAR,nAG,ag+(nAG*(ar-1)));hold on
        x = [];
        for i = 1:size(Rtime_meanDir{ar,ag},3)
        x = vertcat(x,Rtime_meanDir{ar,ag}(:,:,i));
        end
        y = repmat(dirs,nBins,1);
        uY = unique(y);
        [coeff,score] = pca(zscore(x));
        clrs = hsv(length(uY));
        for c = 1:length(uY)
        idx = y==uY(c);
        plot3(score(idx,2),score(idx,3),score(idx,4),'Color',clrs(c,:))
        end
        xlabel('PC2');ylabel('PC3');zlabel('PC4')
        axis square
        box on
        clear score
        title([areas{ar} ' P' num2str(ageGroups{ag}(1)) ' - ' num2str(ageGroups{ag}(2))])
    end
end


figure('Position',[100 100 1000 600])
for ar = 1:nAR
    for ag = 1:nAG
        subplot(nAR,nAG,ag+(nAG*(ar-1)));hold on
        x = Rtime_meanDir{ar,ag};
        for b = 1:nBins
            rB = zscore(x(:,:,b));
            pcorrDir(b) = partialcorr(dir.diss',pdist(rB,'spearman')',ori.diss','Type','Spearman');
            pcorrOri(b) = partialcorr(ori.diss',pdist(rB,'spearman')',dir.diss','Type','Spearman');
        end
        plot(bins,pcorrDir)
        plot(bins,pcorrOri)
        legend({'dir pcorr','ori pcorr'})
        xlabel('time (s)')
        ylabel('partial correlation')
        axis square
        box on
        clear score
        title([areas{ar} ' P' num2str(ageGroups{ag}(1)) ' - ' num2str(ageGroups{ag}(2))])
    end
end

figure('Position',[100 100 1000 600])
for ar = 1:nAR
    for ag = 1:nAG
        subplot(nAR,nAG,ag+(nAG*(ar-1)));hold on
        for b = 1:nBins
            rB = Rtime{ar,ag}(:,:,b);
            accDir(b) = lda_bn(rB,Cdir);
            accOri(b) = lda_bn(rB,Cori);
        end
        plot(bins,accDir)
        plot(bins,accOri)
        legend({'dir','ori'})
        xlabel('time (s)')
        ylabel('LDA decoding accuracy (5-fold cv)')
        axis square
        box on
        clear score
        title([areas{ar} ' P' num2str(ageGroups{ag}(1)) ' - ' num2str(ageGroups{ag}(2))])
    end
end



figure('Position',[100 100 1000 600]);
for ar = 1:nAR
    for ag = 1:nAG

        x = Rtime{ar,ag};
        y = dirs;
        for b = 1:size(x,3)
        for u = 1:size(x,2)
        
            fr = x(:,u,b);
            fr = reshape(fr,nReps,nConds);
            dirDat = getDirTuning(fr,y,0);
            oriDat = getOriTuning(fr,y,0);
        
            dsi(u,b) = dirDat.DSI;
            ldr(u,b) = dirDat.Ldir;
            osi(u,b) = oriDat.OSI;
            lor(u,b) = oriDat.Lori;
        
        end
        end
        idx = bins>=0&bins<=1;

        subplot(nAR,nAG,ag+(nAG*(ar-1)));hold on
        plot(bins(idx),mean(ldr(:,idx),'omitnan'))
        plot(bins(idx),mean(lor(:,idx),'omitnan'))
        legend({'Ldir','Lori'})
        axis square
        box on
        ylim([0 1])
        title([areas{ar} ' P' num2str(ageGroups{ag}(1)) ' - ' num2str(ageGroups{ag}(2))])
    end
end

figure('Position',[100 100 1000 600]);
for ar = 1:nAR
    for ag = 1:nAG

        x = Rtime{ar,ag};
        y = dirs;
        for b = 1:size(x,3)
        for u = 1:size(x,2)
        
            fr = x(:,u,b);
            fr = reshape(fr,nReps,nConds);
            dirDat = getDirTuning(fr,y,0);
            oriDat = getOriTuning(fr,y,0);
        
            dsi(u,b) = dirDat.DSI;
            ldr(u,b) = dirDat.Ldir;
            osi(u,b) = oriDat.OSI;
            lor(u,b) = oriDat.Lori;
        
        end
        end
        idx = bins>=0&bins<=1;

        subplot(nAR,nAG,ag+(nAG*(ar-1)));hold on
        plot(bins(idx),mean(dsi(:,idx),'omitnan'))
        plot(bins(idx),mean(osi(:,idx),'omitnan'))
        legend({'DSI','OSI'})
        axis square
        box on
        ylim([0 1])
        title([areas{ar} ' P' num2str(ageGroups{ag}(1)) ' - ' num2str(ageGroups{ag}(2))])
    end
end

% h=figure;
% colormap gray
% ar = 2;ag = 2;
% A = squareform(pdist(zscore(Rtime_meanDir{ar,ag}(:,:,1)),'spearman'));
% imagesc(A)
% ax = gca;
% ax.NextPlot = 'replaceChildren';
% frames = 61;
% M(frames) = struct('cdata',[],'colormap',[]);
% h.Visible = 'off';
% for f = 1:frames
% B = squareform(pdist(zscore(Rtime_meanDir{ar,ag}(:,:,f)),'spearman'));
% imagesc(B)
% drawnow
% M(f) = getframe;
% end
% h.Visible = 'on';
% movie(M,1,6)
