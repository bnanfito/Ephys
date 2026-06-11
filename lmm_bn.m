% linear mixed effects model - Brandon Nanfito
clear
% close all

% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
dataFold = 'Y:\Brandon\data';
proj = {'Train_BiDir','Train_V1Cool','Train_Static'};
projStyl = {'-','-','-'};
projClr = {'r','g','b'};
anaMode = 'MU';
area = {'V1','PSS'};
exclude = {'febk7','febm6'};

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
        groupId{u} = [num2str(d{p}.manipulation(u)) d{p}.area{u}];
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
    for a = 1:length(area)
%         trainData = d{p}(strcmp(d{p}.area,area{a}),[3,20,22,26:28]);
        trainData = d{p}(strcmp(d{p}.area,area{a}),{'animal','age','area','manipulation','group','ldr','lor'});
        mdl{p,a} = fitlme(trainData,'ldr ~ 1+ manipulation + (1+manipulation|animal)','FitMethod','ML');
    end
end

%% Plot

% cdf plot of tuning metrics
% ldr
figure;hold on
cdfPlt = [];
pltCount = 0;
for p = 1:length(proj)
    pltCount = pltCount+1;
    p1 = cdfplot( d{p}.ldr(strcmp(d{p}.area,'PSS')&d{p}.manipulation==0) );
    p1.Color = projClr{p};
    p1.LineStyle = '--';
    cdfPlt(pltCount) = p1;
    pltLbl{pltCount} = [proj{p}(7:end) ' before training'];
    pltCount = pltCount+1;
    p2 = cdfplot( d{p}.ldr(strcmp(d{p}.area,'PSS')&d{p}.manipulation==1) );
    p2.Color = projClr{p};
    p2.LineStyle = '-';
    cdfPlt(pltCount) = p2;
    pltLbl{pltCount} = [proj{p}(7:end) ' after training'];
end
legend(cdfPlt,pltLbl,'Location','best')

%lor
figure;hold on
cdfPlt = [];
pltCount = 0;
for p = 1:length(proj)
    pltCount = pltCount+1;
    p1 = cdfplot( d{p}.lor(strcmp(d{p}.area,'PSS')&d{p}.manipulation==0) );
    p1.Color = projClr{p};
    p1.LineStyle = '--';
    cdfPlt(pltCount) = p1;
    pltLbl{pltCount} = [proj{p}(7:end) ' before training'];
    pltCount = pltCount+1;
    p2 = cdfplot( d{p}.lor(strcmp(d{p}.area,'PSS')&d{p}.manipulation==1) );
    p2.Color = projClr{p};
    p2.LineStyle = '-';
    cdfPlt(pltCount) = p2;
    pltLbl{pltCount} = [proj{p}(7:end) ' after training'];
end
legend(cdfPlt,pltLbl,'Location','best')

figure
tiledlayout(length(proj),length(area))
for p = 1:length(proj)
    animals = unique(d{p}.animal);
    nAnim = length(animals);
    animClrs = prism(nAnim);
    for ar = 1:length(area)
        nexttile;hold on
        idxArea = strcmp(d{p}.area,area{ar});
        for anim = 1:nAnim
            idxAnim = strcmp(d{p}.animal,animals{anim});
            idxManip = d{p}.manipulation==1;
            idx1 = idxAnim & idxArea & ~idxManip;
            idx2 = idxAnim & idxArea &  idxManip;
            x1 = repmat(anim,sum(idx1),1)+randi([-10 10],sum(idx1),1)/40;
            y1 = d{p}.ldr(idx1);
            x2 = repmat(anim+nAnim+1,sum(idx2),1)+randi([-10 10],sum(idx2),1)/40;
            y2 = d{p}.ldr(idx2);
            plot(x1,y1,'o','Color',animClrs(anim,:))
            plot(x2,y2,'o','Color',animClrs(anim,:))
            xticks([1:nAnim,nAnim+2:nAnim+1+nAnim])
            xticklabels([animals;animals])
        end
        plot([0.5 nAnim+0.5],repmat(mdl{p,ar}.Coefficients.Estimate(1),1,2),'k-','LineWidth',2)
        plot([0.5 nAnim+0.5]+nAnim+1,repmat(mdl{p,ar}.Coefficients.Estimate(1)+mdl{p,ar}.Coefficients.Estimate(2),1,2),'k-','LineWidth',2)
        title([proj{p}(7:end) '; ' area{ar} '; Effect:' num2str(mdl{p,ar}.Coefficients.Estimate(2))])
    end
end



