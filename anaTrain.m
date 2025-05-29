%anaTrain

function [projectTbl,stats,data] = anaTrain(proj)

    anaMode = 'MU';
    sve = 1;
    thresh = 'threshold4';
    visTest = 'ranksum';

    dataFold = 'Y:\Brandon\data';
%     dataFold = '/Volumes/NielsenHome2/Brandon/data';
%     dataFold = fullfile(dataFold,'dataSets','training',proj,anaMode,thresh);

%% load project table

    projectTbl=getProjectFiles(proj,1,'age','recSite','penNr','priorMFlag','priorDescr','duringMFlag','manipDescr','manipDetail');
%     projectTbl = projectTbl(~(strcmp(projectTbl.experimentId,'febm8')|strcmp(projectTbl.experimentId,'febn9')),:);

%     load(fullfile(dataFold,[proj '_' anaMode 'dataSet.mat']),'projectTbl')
%     projectTbl = projectTbl(:,~strcmp(projectTbl.Properties.VariableNames,'sumStats'));

%% Generate SumStats

    for e = 1:height(projectTbl)
        animal = projectTbl.experimentId{e};
        unit = projectTbl.unitNr{e};
        expt = projectTbl.experimentNr{e};
        probe = projectTbl.probeId(e);
        exptName = projectTbl.fileBase{e};
        disp(['generating sumStats for ' exptName])
        [sumStats{e,1}] = anaOri(animal,unit,expt,probe,anaMode,dataFold,0,0);
    end
    projectTbl.sumStats = sumStats;

%% Organize Data
    
    v1bf = vertcat(projectTbl.sumStats{projectTbl.priorMFlag == 0 & strcmp(projectTbl.recSite,'V1')});
%     v1bf = v1bf(screenUnits(v1bf,anaMode),:);
    v1af = vertcat(projectTbl.sumStats{projectTbl.priorMFlag == 1 & strcmp(projectTbl.recSite,'V1')});
%     v1af = v1af(screenUnits(v1af,anaMode),:);
    
    pssbf = vertcat(projectTbl.sumStats{projectTbl.priorMFlag == 0 & strcmp(projectTbl.recSite,'PSS')});
%     pssbf = pssbf(screenUnits(pssbf,anaMode),:);
    pssaf = vertcat(projectTbl.sumStats{projectTbl.priorMFlag == 1 & strcmp(projectTbl.recSite,'PSS')});
%     pssaf = pssaf(screenUnits(pssaf,anaMode),:);

    data.v1bf = v1bf;
    data.v1af = v1af;
    data.pssbf = pssbf;
    data.pssaf = pssaf;

%% Plot

    animals = unique(projectTbl.experimentId);
    aniMarks = {'o','square','diamond','pentagram','^','v','<','>'};
    metrics = {'rPref','dsi','ldr','osi','lor','latency'};
    
    % plot metric cdf across animals
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
            tbl = tbl(screenUnits(tbl,anaMode,visTest),:);
            dist = tbl{:,metrics{m}};
            dist = dist(~isnan(dist));
            N(i) = length(dist);
            DIST{m,i} = dist;

            cdf = cdfplot(dist);
            cdf.Color = clr;
            cdf.LineStyle = linStyl;
            xlabel(metrics{m})
            ylabel('percentile')
            title([])
    
            for a = 1:length(animals)
                
                tmpAniIdx = contains(tbl.exptName,animals{a});
                dist = tbl{tmpAniIdx,metrics{m}};
                dist = dist(~isnan(dist));
                aniDist{m,i,a} = dist;
                aniMean(m,i,a) = mean(dist,'omitnan');
                aniSEM(m,i,a) = std(dist,'omitnan')/sqrt(length(dist));
    
            end
        
        end
        legend({['V1 before; n=' num2str(N(1))],['V1 after; n=' num2str(N(2))],...
                ['PSS before; n=' num2str(N(3))],['PSS after; n=' num2str(N(4))]})
        title([proj ' ' anaMode])
        if sve == 1
            saveas(gcf,fullfile(dataFold,[proj '_' anaMode 'cdf_' metrics{m}]))
        end

    end

    % plot per animal metrics
    for m = 1:length(metrics)

        figure; hold on
        for a = 1:length(animals)
    
            xOff = randn(1)*0.06;
            
            subplot(1,2,1);hold on
            legLbls1Target(a) = plot((1:2)+xOff , aniMean(m,1:2,a),['b-' aniMarks{a}],'LineWidth',1,'MarkerSize',10);
            legLbls1{a} = [animals{a} '; n1=' num2str(length(aniDist{m,1,a})) '; n2=' num2str(length(aniDist{m,2,a}))];
            plot(repmat(1:2,2,1)+xOff , aniMean(m,1:2,a) + ([1;-1]*aniSEM(m,1:2,a)) ,'b-','LineWidth',1)
            xticks(1:2);
            xticklabels({'before','after'})   
            xlim([0.5 2.5])
            ylabel(metrics{m})
            title('V1')

    
            subplot(1,2,2);hold on
            legLbls2Target(a) = plot((1:2)+xOff , aniMean(m,3:4,a),['r-' aniMarks{a}],'LineWidth',1,'MarkerSize',10);
            legLbls2{a} = [animals{a} '; n1=' num2str(length(aniDist{m,3,a})) '; n2=' num2str(length(aniDist{m,4,a}))];
            plot(repmat(1:2,2,1)+xOff , aniMean(m,3:4,a) + ([1;-1]*aniSEM(m,3:4,a)) ,'r-','LineWidth',1)
            xticks(1:2);
            xticklabels({'before','after'})
            xlim([0.5 2.5])
            title('PSS')
    
        end
        subplot(1,2,1);hold on
        legend(legLbls1Target,legLbls1)
        subplot(1,2,2);hold on
        legend(legLbls2Target, legLbls2)
        sgtitle([proj ' ' anaMode])
        if sve == 1
            saveas(gcf,fullfile(dataFold,[proj '_' anaMode 'animals_' metrics{m}]))
        end

    end

%% Stats

    count = 0;
    for m = 1:length(metrics)

        count = count+1;
        aniLbl{count,1} = 'all';
        metLbl{count,1} = metrics{m};
        compLbl{count,1} = 'V1 before vs V1 after';
        tst1Lbl{count,1} = 'kstest2';
        [h1(count,1),p1(count,1)] = kstest2(DIST{m,1},DIST{m,2});
        tst2Lbl{count,1} = 'ranksum';
        [p2(count,1),h2(count,1)] = ranksum(DIST{m,1},DIST{m,2});
        
        count = count+1;
        aniLbl{count,1} = 'all';
        metLbl{count,1} = metrics{m};
        compLbl{count,1} = 'PSS before vs PSS after';
        tst1Lbl{count,1} = 'kstest2';
        [h1(count,1),p1(count,1)] = kstest2(DIST{m,3},DIST{m,4});
        tst2Lbl{count,1} = 'ranksum';
        [p2(count,1),h2(count,1)] = ranksum(DIST{m,3},DIST{m,4});
        

        for a = 1:length(animals)

            count = count+1;
            aniLbl{count,1} = animals{a};
            metLbl{count,1} = metrics{m};
            compLbl{count,1} = 'V1 before vs V1 after';
            tst1Lbl{count,1} = 'kstest2';
            tst2Lbl{count,1} = 'ranksum';
            if ~isempty(aniDist{m,1,a}) && ~isempty(aniDist{m,2,a})
                [h1(count,1),p1(count,1)] = kstest2(aniDist{m,1,a},aniDist{m,2,a});
                [p2(count,1),h2(count,1)] = ranksum(aniDist{m,1,a},aniDist{m,2,a});
            else
                h1(count,1) = false; p1(count,1) = nan;
                h2(count,1) = false; p2(count,1) = nan;
            end
            
            count = count+1;
            aniLbl{count,1} = animals{a};
            metLbl{count,1} = metrics{m};
            compLbl{count,1} = 'PSS before vs PSS after';
            tst1Lbl{count,1} = 'kstest2';
            tst2Lbl{count,1} = 'ranksum';
            if ~isempty(aniDist{m,3,a}) && ~isempty(aniDist{m,4,a})
                [h1(count,1),p1(count,1)] = kstest2(aniDist{m,3,a},aniDist{m,4,a});
                [p2(count,1),h2(count,1)] = ranksum(aniDist{m,3,a},aniDist{m,4,a});
            else
                h1(count,1) = false; p1(count,1) = nan;
                h2(count,1) = false; p2(count,1) = nan;
            end

        end

    end
    varNames = {'animal','metric','comparison','test1','h1','p1','test2','h2','p2'};
    stats = table(aniLbl,metLbl,compLbl,tst1Lbl,h1,p1,tst2Lbl,h2,p2,'VariableNames',varNames);

    if sve == 1
        save(fullfile(dataFold,[proj '_' anaMode 'dataSet.mat']),'projectTbl','data','stats')
    end

end