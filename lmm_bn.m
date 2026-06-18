% linear mixed effects model - Brandon Nanfito
clear
close all

% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
dataFold = 'Y:\Brandon\data';
proj = {'Train_BiDir','Train_V1Cool','Train_Static'};
projStyl = {'-','-','-'};
projShap = {'o','v','square'};
projClr = {'r','g','b'};
anaMode = 'MU';
area = {'V1','PSS'};
exclude = {};

%% Organize data

for p = 1:length(proj)
    load(fullfile(dataFold,'dataSets','training',proj{p},anaMode,[proj{p} '_' anaMode 'dataSet.mat']),'-mat','data','projectTbl');
    d{p} = data;
    manipFlag = [zeros(height(d{p}.v1bf),1);ones(height(d{p}.v1af),1);...
                 zeros(height(d{p}.pssbf),1);ones(height(d{p}.pssaf),1)];
    d{p} = vertcat(d{p}.v1bf,d{p}.v1af,d{p}.pssbf,d{p}.pssaf);
    d{p}.manipulation = manipFlag;
    aniId = cell(height(d{p}),1);
    uAge = nan(height(d{p}),1);
    groupId = cell(height(d{p}),1);
    for u = 1:height(d{p})
        aniId{u} = d{p}.exptName{u}(1:5);
        uAge(u) = unique( projectTbl.age(strcmp(projectTbl.experimentId,aniId{u})) );
        groupId{u} = [proj{p} num2str(d{p}.manipulation(u)) d{p}.area{u}];
    end
    d{p}.animal = aniId;
    d{p}.age = uAge;
    d{p}.group = groupId;
    d{p} = d{p}(~ismember(d{p}.animal,exclude),:);
    d{p} = d{p}(screenUnits(d{p},'MU'),:);
    clear data projectTbl aniId uAge manipFlag
end

%% LMM

for p = 1:length(proj)
    for ar = 1:length(area)
%         trainData = d{p}(strcmp(d{p}.area,area{a}),[3,20,22,26:28]);
        trainData = d{p}(strcmp(d{p}.area,area{ar}),{'animal','age','area','manipulation','group','ldr','lor'});
        
        mdl_ldr{p,ar} = fitlme(trainData,'ldr ~ 1+ manipulation + (1+manipulation|animal)','FitMethod','ML');
        ldrInt(p,ar) = mdl_ldr{p,ar}.Coefficients.Estimate(1);
        ldrSlope(p,ar) = mdl_ldr{p,ar}.Coefficients.Estimate(2);
        [~,~,anim_ldr]=randomEffects(mdl_ldr{p,ar});
        anim_ldrName{p,ar} = anim_ldr.Level(strcmp(anim_ldr.Name,'(Intercept)'));
        anim_ldrInt{p,ar} = anim_ldr.Estimate(strcmp(anim_ldr.Name,'(Intercept)'))+ldrInt(p,ar);
        anim_ldrInt_sem{p,ar} = anim_ldr.SEPred(strcmp(anim_ldr.Name,'(Intercept)'));
        anim_ldrSlope{p,ar} = anim_ldr.Estimate(strcmp(anim_ldr.Name,'manipulation'))+ldrSlope(p,ar);
        anim_ldrSlope_sem{p,ar} = anim_ldr.SEPred(strcmp(anim_ldr.Name,'manipulation'));
        anim_ldrPost{p,ar} = anim_ldrInt{p,ar}+anim_ldrSlope{p,ar};

        mdl_lor{p,ar} = fitlme(trainData,'lor ~ 1+ manipulation + (1+manipulation|animal)','FitMethod','ML');
        lorInt(p,ar) = mdl_lor{p,ar}.Coefficients.Estimate(1);
        lorSlope(p,ar) = mdl_lor{p,ar}.Coefficients.Estimate(2);
        [~,~,anim_lor]=randomEffects(mdl_lor{p,ar});
        anim_lorName{p,ar} = anim_lor.Level(strcmp(anim_lor.Name,'(Intercept)'));
        anim_lorInt{p,ar} = anim_lor.Estimate(strcmp(anim_lor.Name,'(Intercept)'))+lorInt(p,ar);
        anim_lorInt_sem{p,ar} = anim_lor.SEPred(strcmp(anim_lor.Name,'(Intercept)'));
        anim_lorSlope{p,ar} = anim_lor.Estimate(strcmp(anim_lor.Name,'manipulation'))+lorSlope(p,ar);
        anim_lorSlope_sem{p,ar} = anim_lor.SEPred(strcmp(anim_lor.Name,'manipulation'));
        anim_lorPost{p,ar} = anim_lorInt{p,ar}+anim_lorSlope{p,ar};

    end
end

%% Plot

% cdf plot of tuning metrics
% ldr
figure('Position',[0 100 1200 800])
tiledlayout(1,length(area))
for ar = 1:length(area)
    nexttile; hold on
    cdfPlt = [];
    pltCount = 0;
    for p = 1:length(proj)
        pltCount = pltCount+1;
        p1 = cdfplot( d{p}.ldr(strcmp(d{p}.area,area{ar})&d{p}.manipulation==0) );
        p1.Color = projClr{p};
        p1.LineStyle = '--';
        cdfPlt(pltCount) = p1;
        pltLbl{pltCount} = [proj{p}(7:end) ' before training'];
        pltCount = pltCount+1;
        p2 = cdfplot( d{p}.ldr(strcmp(d{p}.area,area{ar})&d{p}.manipulation==1) );
        p2.Color = projClr{p};
        p2.LineStyle = '-';
        cdfPlt(pltCount) = p2;
        pltLbl{pltCount} = [proj{p}(7:end) ' after training'];
    end
    xlabel('Ldir')
    ylabel('Proportion')
    legend(cdfPlt,pltLbl,'Location','best')
    axis square; box on
    title(area{ar})
end

figure('Position',[0 100 1200 800])
% tiledlayout(length(proj),length(area))
for p = 1:length(proj)
    animals = unique(d{p}.animal);
    nAnim = length(animals);
    animClrs = parula(nAnim);
    for ar = 1:length(area)
%         nexttile;hold on
        subplot(length(area),length(proj), p+(length(proj)*(ar-1)) );hold on
        idxArea = strcmp(d{p}.area,area{ar});
        for an = 1:nAnim
            idxAnim = strcmp(d{p}.animal,animals{an});
            idxManip = d{p}.manipulation==1;
            idx1 = idxAnim & idxArea & ~idxManip;
            idx2 = idxAnim & idxArea &  idxManip;
            x1 = repmat(an,sum(idx1),1)+randi([-10 10],sum(idx1),1)/40;
            y1 = d{p}.ldr(idx1);
            x2 = repmat(an+nAnim+1,sum(idx2),1)+randi([-10 10],sum(idx2),1)/40;
            y2 = d{p}.ldr(idx2);
            scatter(x1,y1,'o','MarkerEdgeColor','none','MarkerFaceColor',animClrs(an,:),'MarkerFaceAlpha',0.35)
            scatter(x2,y2,'o','MarkerEdgeColor','none','MarkerFaceColor',animClrs(an,:),'MarkerFaceAlpha',0.35)
            xticks([1:nAnim,nAnim+2:nAnim+1+nAnim])
            xticklabels([animals;animals])
        end

        plot([0.5 nAnim+0.5],repmat(ldrInt(p,ar),1,2),'k-','LineWidth',2)
        plot([0.5 nAnim+0.5]+nAnim+1,repmat(ldrInt(p,ar)+ldrSlope(p,ar),1,2),'k-','LineWidth',2)
        
        for an = 1:nAnim
            anim_idx = find(strcmp(animals{an},anim_ldrName{p,ar}));
            if ~isempty(anim_idx)
                plot(an+[-0.5 0.5],repmat(anim_ldrInt{p,ar}(anim_idx),1,2),'r-','LineWidth',2)
                plot(an+nAnim+1+[-0.5 0.5],repmat(anim_ldrPost{p,ar}(anim_idx),1,2),'r-','LineWidth',2)
            end
        end

        title([proj{p}(7:end) '; ' area{ar} '; Effect: ' num2str(ldrSlope(p,ar))])
        ylabel('Ldir')
%         ylim([0 1])
        xlim([0 (2*nAnim)+2])
        box on
    end
end


%lor
figure('Position',[0 100 1200 800])
tiledlayout(1,length(area))
for ar = 1:length(area)
    nexttile; hold on
    cdfPlt = [];
    pltCount = 0;
    for p = 1:length(proj)
        pltCount = pltCount+1;
        p1 = cdfplot( d{p}.lor(strcmp(d{p}.area,area{ar})&d{p}.manipulation==0) );
        p1.Color = projClr{p};
        p1.LineStyle = '--';
        cdfPlt(pltCount) = p1;
        pltLbl{pltCount} = [proj{p}(7:end) ' before training'];
        pltCount = pltCount+1;
        p2 = cdfplot( d{p}.lor(strcmp(d{p}.area,area{ar})&d{p}.manipulation==1) );
        p2.Color = projClr{p};
        p2.LineStyle = '-';
        cdfPlt(pltCount) = p2;
        pltLbl{pltCount} = [proj{p}(7:end) ' after training'];
    end
    xlabel('Lori')
    ylabel('Proportion')
    legend(cdfPlt,pltLbl,'Location','best')
    axis square; box on
    title(area{ar})
end

figure('Position',[0 100 1200 800])
% tiledlayout(length(proj),length(area))
for p = 1:length(proj)
    animals = unique(d{p}.animal);
    nAnim = length(animals);
    animClrs = parula(nAnim);
    for ar = 1:length(area)
%         nexttile;hold on
        subplot(length(area),length(proj), p+(length(proj)*(ar-1)) );hold on
        idxArea = strcmp(d{p}.area,area{ar});
        
        for an = 1:nAnim
            idxAnim = strcmp(d{p}.animal,animals{an});
            idxManip = d{p}.manipulation==1;
            idx1 = idxAnim & idxArea & ~idxManip;
            idx2 = idxAnim & idxArea &  idxManip;
            x1 = repmat(an,sum(idx1),1)+randi([-10 10],sum(idx1),1)/40;
            y1 = d{p}.lor(idx1);
            x2 = repmat(an+nAnim+1,sum(idx2),1)+randi([-10 10],sum(idx2),1)/40;
            y2 = d{p}.lor(idx2);
            scatter(x1,y1,'o','MarkerEdgeColor','none','MarkerFaceColor',animClrs(an,:),'MarkerFaceAlpha',0.35)
            scatter(x2,y2,'o','MarkerEdgeColor','none','MarkerFaceColor',animClrs(an,:),'MarkerFaceAlpha',0.35)
            xticks([1:nAnim,nAnim+2:nAnim+1+nAnim])
            xticklabels([animals;animals])
        end
        
        plot([0.5 nAnim+0.5],repmat(lorInt(p,ar),1,2),'k-','LineWidth',2)
        plot([0.5 nAnim+0.5]+nAnim+1,repmat(lorInt(p,ar)+lorSlope(p,ar),1,2),'k-','LineWidth',2)
        
        for an = 1:nAnim
            anim_idx = find(strcmp(animals{an},anim_lorName{p,ar}));
            if ~isempty(anim_idx)
                plot(an+[-0.5 0.5],repmat(anim_lorInt{p,ar}(anim_idx),1,2),'r-','LineWidth',2)
                plot(an+nAnim+1+[-0.5 0.5],repmat(anim_lorPost{p,ar}(anim_idx),1,2),'r-','LineWidth',2)
            end
        end
       
        title([proj{p}(7:end) '; ' area{ar} '; Effect: ' num2str(lorSlope(p,ar))])
        ylabel('Lori')
%         ylim([0 1])
        xlim([0 (2*nAnim)+2])
        box on
    end
end


figure('Position',[0 100 700 1160])
for p = 1:length(proj)
    for ar = 1:length(area)
        subplot( length(area)+1,2,1+(length(area)*(ar-1)) ); hold on
        scatter(anim_ldrInt{p,ar},anim_ldrSlope{p,ar},projShap{p},'filled')
        yline(0,':')
        xlabel('Pretraining Ldir')
        ylabel('Change in Ldir')
        box on; axis square
        title(area{ar})

        subplot( length(area)+1,2,2+(length(area)*(ar-1)) ); hold on
        scatter(anim_lorInt{p,ar},anim_lorSlope{p,ar},projShap{p},'filled')
        yline(0,':')
        xlabel('Pretraining Lori')
        ylabel('Change in Lori')
        box on; axis square
        title(area{ar})
    end
end

% figure;hold on
for p = 1:length(proj)
    subplot(3,2,5); hold on
    if p==1
        xline(0,'k:')
        yline(0,'k:')
    end
    errorbar(anim_ldrSlope{p,1},anim_ldrSlope{p,2},...
        anim_ldrSlope_sem{p,2},anim_ldrSlope_sem{p,2},...
        anim_ldrSlope_sem{p,1},anim_ldrSlope_sem{p,1},projShap{p},'LineStyle','none','MarkerFaceColor','auto')
    xlabel(area{1})
    ylabel(area{2})
    title('Change in Ldir')
    box on; axis square

    subplot(3,2,6); hold on
    if p==1
        xline(0,'k:')
        yline(0,'k:')
    end
    errorbar(anim_lorSlope{p,1},anim_lorSlope{p,2},...
        anim_lorSlope_sem{p,2},anim_lorSlope_sem{p,2},...
        anim_lorSlope_sem{p,1},anim_lorSlope_sem{p,1},projShap{p},'LineStyle','none','MarkerFaceColor','auto')
    xlabel(area{1})
    ylabel(area{2})
    title('Change in Lori')
    box on; axis square
end

