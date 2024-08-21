%anaTrain_BiDir

clear all
close all

% load('/Volumes/Lab drive/Brandon/data/dataSets/training/Train_BiDir/projectTbl.mat')

projectTbl=getProjectFiles('Train_BiDir',1,'age','recSite','priorMFlag','priorDescr','duringMFlag','manipDescr','manipDetail');
dataFold = 'Y:\Brandon\data';

for e = 1:height(projectTbl)
    animal = projectTbl.experimentId{e};
    unit = projectTbl.unitNr{e};
    expt = projectTbl.experimentNr{e};
    probe = projectTbl.probeId(e);
    exptName = projectTbl.fileBase{e};
    exptDir = fullfile(dataFold,'Ephys',animal,exptName);
    disp(['generating sumStats for ' exptName])
    [sumStats{e,1}] = anaOri(animal,unit,expt,probe,'MU',dataFold,0,0);
end
projectTbl.sumStats = sumStats;

v1bf = vertcat(projectTbl.sumStats{projectTbl.priorMFlag == 0 & strcmp(projectTbl.recSite,'V1')});
v1af = vertcat(projectTbl.sumStats{projectTbl.priorMFlag == 1 & strcmp(projectTbl.recSite,'V1')});

pssbf = vertcat(projectTbl.sumStats{projectTbl.priorMFlag == 0 & strcmp(projectTbl.recSite,'PSS')});
pssaf = vertcat(projectTbl.sumStats{projectTbl.priorMFlag == 1 & strcmp(projectTbl.recSite,'PSS')});

animals = unique(projectTbl.experimentId);
aniMarks = {'o','+','*','x','square','diamond','^','v','<','>'};
metrics = {'rPref','dsi','ldr'};

for m = 1:length(metrics)
    figure; hold on
    for i = 1:4
    
        if i == 1
            area = 'V1';
            tbl = v1bf;
            clr = 'b';
            linStyl = '--';
        elseif i == 2
            area = 'V1';
            tbl = v1af;
            clr = 'b';
            linStyl = '-';
        elseif i == 3
            area = 'PSS';
            tbl = pssbf;
            clr = 'r';
            linStyl = '--';
        elseif i == 4
            area = 'PSS';
            tbl = pssaf;
            clr = 'r';
            linStyl = '-';
        end
    
    %     goodIdx = tbl.goodUnit;
        for u = 1:height(tbl)
    %         pVis(u) = anova1([tbl.response{u} tbl.rBlank{u}],[],"off");
            baseFR = tbl.fr(u).base(:,1:end-1);
            stimFR = tbl.fr(u).stim(:,1:end-1);
            pVis(u) = ranksum(baseFR(:),stimFR(:));
            clear baseFR stimFR
        end
        isVis = pVis<0.01;
        isAct = tbl.rPref>2;
        goodIdx = isVis' & isAct;
        tbl = tbl(goodIdx,:);
        clear pVis
    
        dist = tbl{:,metrics{m}};
        cdf = cdfplot(dist);
        cdf.Color = clr;
        cdf.LineStyle = linStyl;
        xlabel(metrics{m})
        ylabel('percentile')
        title([])

        for a = 1:length(animals)
            
            tmpAniIdx = contains(tbl.exptName,animals{a});
            dist = tbl{tmpAniIdx,metrics{m}};
            aniMean(m,i,a) = mean(dist);
            aniSEM(m,i,a) = std(dist)/sqrt(length(dist));

        end
    
    end
end

for m = 1:length(metrics)
    figure; hold on
    for a = 1:length(animals)

        xOff = randn(1)*0.06;
        
        subplot(1,2,1);hold on
        plot((1:2)+xOff , aniMean(m,1:2,a),['b-' aniMarks{a}],'LineWidth',1,'MarkerSize',10)
        plot(repmat(1:2,2,1)+xOff , aniMean(m,1:2,a) + ([1;-1]*aniSEM(m,1:2,a)) ,'b-','LineWidth',1)
        xticks(1:2);
        xticklabels({'before','after'})   
        xlim([0.5 2.5])
        ylabel(metrics{m})

        subplot(1,2,2);hold on
        plot((1:2)+xOff , aniMean(m,3:4,a),['r-' aniMarks{a}],'LineWidth',1,'MarkerSize',10)
        plot(repmat(1:2,2,1)+xOff , aniMean(m,3:4,a) + ([1;-1]*aniSEM(m,3:4,a)) ,'r-','LineWidth',1)
        xticks(1:2);
        xticklabels({'before','after'})
        xlim([0.5 2.5])

    end
end



