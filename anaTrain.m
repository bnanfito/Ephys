%anaTrain

function [projectTbl,stats] = anaTrain(proj)

%     dataFold = ['/Volumes/Lab drive/Brandon/data/dataSets/training/' proj];
%     load(fullfile(dataFold,'projectTbl.mat'),'projectTbl')

    dataFold = 'Y:\Brandon\data';
    projectTbl=getProjectFiles(proj,1,'age','recSite','priorMFlag','priorDescr','duringMFlag','manipDescr','manipDetail');
    
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
    metrics = {'rPref','dsi','ldr','latency'};

    %% plot
    
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
                aniDist{m,i,a} = dist;
                aniMean(m,i,a) = mean(dist,'omitnan');
                aniSEM(m,i,a) = std(dist,'omitnan')/sqrt(length(dist));
    
            end
        
        end
        legend({'V1 before','V1 after','PSS before','PSS after'})
    end
    
    for m = 1:length(metrics)
        figure; hold on
        for a = 1:length(animals)
    
            xOff = randn(1)*0.06;
            
            subplot(1,2,1);hold on
            legLbls1Target(a) = plot((1:2)+xOff , aniMean(m,1:2,a),['b-' aniMarks{a}],'LineWidth',1,'MarkerSize',10);
            legLbls1{a} = animals{a};
            plot(repmat(1:2,2,1)+xOff , aniMean(m,1:2,a) + ([1;-1]*aniSEM(m,1:2,a)) ,'b-','LineWidth',1)
            xticks(1:2);
            xticklabels({'before','after'})   
            xlim([0.5 2.5])
            ylabel(metrics{m})

    
            subplot(1,2,2);hold on
            legLbls2Target(a) = plot((1:2)+xOff , aniMean(m,3:4,a),['r-' aniMarks{a}],'LineWidth',1,'MarkerSize',10);
            legLbls2{a} = animals{a};
            plot(repmat(1:2,2,1)+xOff , aniMean(m,3:4,a) + ([1;-1]*aniSEM(m,3:4,a)) ,'r-','LineWidth',1)
            xticks(1:2);
            xticklabels({'before','after'})
            xlim([0.5 2.5])
    
        end
        subplot(1,2,1);hold on
        legend(legLbls1Target,legLbls1)
        subplot(1,2,2);hold on
        legend(legLbls2Target, legLbls2)
    end

    %% stats

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


end