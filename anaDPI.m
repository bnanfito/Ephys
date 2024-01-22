 %anaDPI
%analyze the hour by hour changes in direction preference index of units
%during training

clear all
close all

if ispc
%     dataFold = 'D:\data';
    dataFold = 'D:\OneDrive - Johns Hopkins\Documents\data';
elseif ismac
    dataFold = '/Volumes/Lab drive/Brandon/data';
end
physDir = fullfile(dataFold,'Ephys');
figDir = fullfile(dataFold,'Figures');
sumDir = fullfile(dataFold,'SummaryStats');
anaMode = 'MU';

animals = {'FEAO4','FEAN6','FEAS6','FEAT1','FEAQ5'};
% animals = {'febj4'};
% animals = {'febh2','febh3','febj3'};
nAnimals = length(animals);
figure;hold on
for a = 1:nAnimals
    animal = animals{a};
    
    if strcmp(animal,'FEAO4')
    
        blocks{1,a} = 'FEAO4_u001_000';
        blocks{2,a} = 'FEAO4_u001_001';
        blocks{3,a} = 'FEAO4_u001_002';
        blocks{4,a} = 'FEAO4_u001_003';
        blocks{5,a} = 'FEAO4_u001_004';
        blocks{6,a} = 'FEAO4_u001_005';
        blocks{7,a} = 'FEAO4_u001_006';
        blocks{8,a} = 'FEAO4_u001_007';
        blocks{9,a} = 'FEAO4_u001_008';
        blocks{10,a} = 'FEAO4_u001_009';
        blocks{11,a} = 'FEAO4_u001_010';
        blocks{12,a} = 'FEAO4_u001_011';
        blocks{13,a} = 'FEAO4_u001_012';
        blocks{14,a} = 'FEAO4_u001_013';
        blocks{15,a} = 'FEAO4_u001_015';
        blocks{16,a} = 'FEAO4_u001_016';
    
        flagBlocks = [];
    
    elseif strcmp(animal,'FEAN6')
    
        blocks{1,a} = 'FEAN6_u001_000';
        blocks{2,a} = 'FEAN6_u001_001';
        blocks{3,a} = 'FEAN6_u001_002';
        blocks{4,a} = 'FEAN6_u001_003';
        blocks{5,a} = 'FEAN6_u001_004';
        blocks{6,a} = 'FEAN6_u001_005';
        blocks{7,a} = 'FEAN6_u001_006';
        blocks{8,a} = 'FEAN6_u001_007';
        blocks{9,a} = 'FEAN6_u001_008';
        blocks{10,a} = 'FEAN6_u001_009';
        blocks{11,a} = 'FEAN6_u001_010';
        blocks{12,a} = 'FEAN6_u001_011';
        blocks{13,a} = 'FEAN6_u001_012';
        blocks{14,a} = 'FEAN6_u001_013';
        blocks{15,a} = 'FEAN6_u001_014';
        blocks{16,a} = 'FEAN6_u001_015';
    
        flagBlocks = [];
    
    elseif strcmp(animal,'FEAS6')
    
        blocks{1,a} = 'FEAS6_u001_000';
        blocks{2,a} = 'FEAS6_u001_001';
        blocks{3,a} = 'FEAS6_u001_002';
        blocks{4,a} = 'FEAS6_u001_003';
        blocks{5,a} = 'FEAS6_u001_004';
        blocks{6,a} = 'FEAS6_u001_005';
        blocks{7,a} = 'FEAS6_u001_006';
        blocks{8,a} = 'FEAS6_u001_007';
        blocks{9,a} = 'FEAS6_u001_008';
        blocks{10,a} = 'FEAS6_u001_009'; %digi inputs missing
        blocks{11,a} = 'FEAS6_u001_010'; %digi inputs missing
        blocks{12,a} = 'FEAS6_u001_011'; %digi inputs missing
        blocks{13,a} = 'FEAS6_u001_012'; %digi inputs missing
        blocks{14,a} = 'FEAS6_u001_013'; %digi inputs missing
        blocks{15,a} = 'FEAS6_u001_014';
        blocks{16,a} = 'FEAS6_u001_015';
    
        flagBlocks = 10:14;
    
    elseif strcmp(animal,'FEAT1')
    
        blocks{1,a} = 'FEAT1_u001_000';
        blocks{2,a} = 'FEAT1_u001_001';
        blocks{3,a} = 'FEAT1_u001_002';
        blocks{4,a} = 'FEAT1_u001_003';
        blocks{5,a} = 'FEAT1_u001_004';
        blocks{6,a} = 'FEAT1_u001_005';
        blocks{7,a} = 'FEAT1_u001_006';
        blocks{8,a} = 'FEAT1_u001_007';
        blocks{9,a} = 'FEAT1_u001_008';
        blocks{10,a} = 'FEAT1_u001_009';
        blocks{11,a} = 'FEAT1_u001_010';
        blocks{12,a} = 'FEAT1_u001_011';
        blocks{13,a} = 'FEAT1_u001_012';
        blocks{14,a} = 'FEAT1_u001_013';
        blocks{15,a} = 'FEAT1_u001_014';
        blocks{16,a} = 'FEAT1_u001_015';
    
        flagBlocks = [];
    
    elseif strcmp(animal,'FEAQ5')
    
        blocks{1,a} = 'FEAQ5_u001_000';
        blocks{2,a} = 'FEAQ5_u001_001';
        blocks{3,a} = 'FEAQ5_u001_002';
        blocks{4,a} = 'FEAQ5_u001_003';
        blocks{5,a} = 'FEAQ5_u001_004';
        blocks{6,a} = 'FEAQ5_u001_005';
        blocks{7,a} = 'FEAQ5_u001_006';
        blocks{8,a} = 'FEAQ5_u001_007';
        blocks{9,a} = 'FEAQ5_u001_008';
        blocks{10,a} = 'FEAQ5_u001_009';
        blocks{11,a} = 'FEAQ5_u001_010';
        blocks{12,a} = 'FEAQ5_u001_011';
        blocks{13,a} = 'FEAQ5_u001_012';
        blocks{14,a} = 'FEAQ5_u001_013';
        blocks{15,a} = 'FEAQ5_u001_014';
        blocks{16,a} = 'FEAQ5_u001_015';
    
        flagBlocks = [];
    
    elseif strcmp(animal,'febh2')

        blocks{1,a} = 'febh2_u000_005'; %electrical noise issues (from pump)
        blocks{2,a} = 'febh2_u000_007'; %electrical noise issues (from pump)
        blocks{3,a} = 'febh2_u000_009';
        blocks{4,a} = 'febh2_u000_011';
        blocks{5,a} = 'febh2_u000_013';
        blocks{6,a} = 'febh2_u000_016';
        blocks{7,a} = 'febh2_u000_018'; %digi inputs missing
        blocks{8,a} = 'febh2_u000_020'; %digi inputs missing
        blocks{9,a} = 'febh2_u000_022'; %digi inputs missing
        blocks{10,a} = 'febh2_u000_024'; %digi inputs missing
        blocks{11,a} = 'febh2_u000_026'; %digi inputs missing
        blocks{12,a} = 'febh2_u000_028';
        blocks{13,a} = 'febh2_u000_030';
        blocks{14,a} = 'febh2_u000_032';
        blocks{15,a} = 'febh2_u000_034';
        blocks{16,a} = 'febh2_u000_036';

        flagBlocks = [];

    elseif strcmp(animal,'febh3')

        blocks{1,a} = 'febh3_u000_006';
        blocks{2,a} = 'febh3_u000_008';
        blocks{3,a} = 'febh3_u000_010';
        blocks{4,a} = 'febh3_u000_014'; %something weird with the trialInfo file
        blocks{5,a} = 'febh3_u000_017';
        blocks{6,a} = 'febh3_u000_019';
        blocks{7,a} = 'febh3_u000_021';
        blocks{8,a} = 'febh3_u000_023';
        blocks{9,a} = 'febh3_u000_025';
        blocks{10,a} = 'febh3_u000_027';
        blocks{11,a} = 'febh3_u000_029';
        blocks{12,a} = 'febh3_u000_031';
        blocks{13,a} = 'febh3_u000_033';
        blocks{14,a} = 'febh3_u000_035';
        blocks{15,a} = 'febh3_u000_037';
        blocks{16,a} = 'febh3_u000_039';

        flagBlocks = [4];

    elseif strcmp(animal,'febj3')

        blocks{1,a} = 'febj3_u000_005';
        blocks{2,a} = 'febj3_u000_007';
        blocks{3,a} = 'febj3_u000_009';
        blocks{4,a} = 'febj3_u000_011';
        blocks{5,a} = 'febj3_u000_013';
        blocks{6,a} = 'febj3_u000_015';
        blocks{7,a} = 'febj3_u000_017';
        blocks{8,a} = 'febj3_u000_019';
        blocks{9,a} = 'febj3_u000_021';
        blocks{10,a} = 'febj3_u000_023';
        blocks{11,a} = 'febj3_u000_025';
        blocks{12,a} = 'febj3_u000_028'; %fewer reps, block split with febj3_u000_027
        blocks{13,a} = 'febj3_u000_030';
        blocks{14,a} = 'febj3_u000_032';
        blocks{15,a} = 'febj3_u000_034';
        blocks{16,a} = 'febj3_u000_036';

        flagBlocks = [];

    elseif strcmp(animal,'febj4')

        blocks{1,a} = 'febj4_u000_008';
        blocks{2,a} = 'febj4_u000_009';
        blocks{3,a} = 'febj4_u000_010';
        blocks{4,a} = 'febj4_u000_011';
        blocks{5,a} = 'febj4_u000_012';
        blocks{6,a} = 'febj4_u000_013';
        blocks{7,a} = 'febj4_u000_014';
        blocks{8,a} = 'febj4_u000_015';
        blocks{9,a} = 'febj4_u000_016';
        blocks{10,a} = 'febj4_u000_017';
        blocks{11,a} = 'febj4_u000_018';
        blocks{12,a} = 'febj4_u000_019';
        blocks{13,a} = 'febj4_u000_020';
        blocks{14,a} = 'febj4_u000_021';
        blocks{15,a} = 'febj4_u000_022';
        blocks{16,a} = 'febj4_u000_023';

        flagBlocks = [];

    end
    nBlocks = size(blocks,1);
    nHr = nBlocks/2;
    
    for b = 1:nBlocks
    
        exptName = blocks{b,a};
        unit = exptName(8:10);
        expt = exptName(12:14);
    
        if ismember(b,flagBlocks)
    
            disp([exptName ' skipped'])
    
        else
    
            disp(exptName)
            load(fullfile(physDir,animal,exptName,[exptName '_id.mat']),'id')
            load(fullfile(physDir,animal,exptName,[exptName '_trialInfo.mat']),'trialInfo')
            for p = 1:length(id.probes)
                sumFile = fullfile(sumDir,animal,exptName,[exptName '_p' num2str(p) '_sumStats' anaMode '.mat']);
                if isfile(sumFile)
                    load(sumFile,'sumStats')
                else
                    [sumStats,spks] = plotUnits(animal,unit,expt,p,anaMode,'ranksum',0.01,0,0,0,dataFold);
                end
                if strcmp(id.probes(p).area,'V1')
                    v1{b,a} = sumStats;
                elseif strcmp(id.probes(p).area,'PSS')
                    pss{b,a} = sumStats;
                end
                clear sumStats
            end
    
        end
    
    end
    clear b

end

animals{nAnimals+1} = 'all Animals';
for b = 1:nBlocks
    v1{b,nAnimals+1} = vertcat(v1{b,1:nAnimals});
    pss{b,nAnimals+1} = vertcat(pss{b,1:nAnimals});
end


for a = 1:length(animals)

    dpiV1{a} = [];
    dpiPSS{a} = [];

    for p = 1:2
    
        countLbl = 0;

        if p ==1 %index for V1
            sp = subplot(nAnimals+1,2,1+(2*(a-1)));hold on
        elseif p == 2 % index for PSS
            sp = subplot(nAnimals+1,2,2+(2*(a-1)));hold on
        end

        for b = 1:2:nBlocks

            hr = ceil(b/2);
            color = (hr/(nHr));

            if p==1
                color = [0 0 color];
                area = 'V1';
                tbl = vertcat(v1{b,a},v1{b+1,a});
            elseif p==2
                color = [color 0 0];
                area = 'PSS';
                tbl = vertcat(pss{b,a},pss{b+1,a});
            end
    
            if isempty(tbl) || sum(tbl.goodUnit)==0
                continue
            else
                countLbl = countLbl+1;
            end

            dist = tbl.dpi(tbl.goodUnit==1);
            if strcmp(area,'V1')
                dpiV1{a} = vertcat(dpiV1{a},[dist repmat(hr,length(dist),1)]);
            elseif strcmp(area,'PSS')
                dpiPSS{a} = vertcat(dpiPSS{a},[dist repmat(hr,length(dist),1)]);
            end
            c = cdfplot(dist);
            c.Color = color;
            c.LineWidth = 2;
            yticks([0 0.5 1])
            ylabel(animals{a})
            xlim([0 1])
            xlabel('DPI')
            lbl{countLbl} = ['hour' num2str(hr)];

            if strcmp(area,'V1')
                v1Violin{a,hr} = violin(dist,0.3,[0 1 1000],1,0);
                medV1(a,hr) = median(dist);
                meanV1(a,hr) = mean(dist);
            elseif strcmp(area,'PSS')
                pssViolin{a,hr} = violin(dist,0.3,[0 1 1000],1,0);
                medPSS(a,hr) = median(dist);
                meanPSS(a,hr) = mean(dist);
            end
    
        end
        legend(lbl)
    
    end

end

figure;hold on
for hr = 1:hr

    subplot(2,1,1);hold on
    x = v1Violin{end,hr}.kdeX;
    y = v1Violin{end,hr}.kdeY;
    plot(hr+([x*0.45 fliplr(x)*-0.45]),[y fliplr(y)],'b','LineWidth',2)
    patch(hr+([x*0.45 fliplr(x)*-0.45]),[y fliplr(y)],'b','EdgeColor','none','FaceAlpha',0.2)
    clear x y
    plot(hr+[0.5,-0.5],repmat(meanV1(end,hr),1,2),'k','LineWidth',2)
    plot(hr+[0.5,-0.5],repmat(medV1(end,hr),1,2),'Color',[0.5 0.5 0.5],'LineWidth',2)
    ylabel('DPI')

    subplot(2,1,2);hold on
    x = pssViolin{end,hr}.kdeX;
    y = pssViolin{end,hr}.kdeY;
    plot(hr+([x*0.45 fliplr(x)*-0.45]),[y fliplr(y)],'r','LineWidth',2)
    patch(hr+([x*0.45 fliplr(x)*-0.45]),[y fliplr(y)],'r','EdgeColor','none','FaceAlpha',0.2)
    clear x y
    plot(hr+[0.5,-0.5],repmat(meanPSS(end,hr),1,2),'k','LineWidth',2)
    plot(hr+[0.5,-0.5],repmat(medPSS(end,hr),1,2),'Color',[0.5 0.5 0.5],'LineWidth',2)
    xlabel('hours of stimulation')
    ylabel('DPI')

end


figure;hold on
hrs = [1 4 8];
shapes = {'o','square','diamond','^','v'};
markerSize = 7;
for a = 1:nAnimals
    subplot(1,2,1);hold on
    y = medV1(a,hrs)';
    plot(repmat(hrs,nAnimals,1)',y,'b-','LineWidth',2,'Marker',shapes{a},'MarkerSize',markerSize)
    xticks(hrs)
    ylim([0 1])
    ylabel('DPI')
    xlabel('hours or stimulation')
    
    subplot(1,2,2);hold on
    y = medPSS(a,hrs)';
    plot(repmat(hrs,nAnimals,1)',y,'r-','LineWidth',2,'Marker',shapes{a},'MarkerSize',markerSize)
    xticks(hrs)
    ylim([0 1])
    xlabel('hours or stimulation')
    sgtitle('median DPI (per animal)')
end



figure;hold on
for a = 1:nAnimals
    subplot(1,2,1);hold on
    y = meanV1(a,hrs)';
    plot(repmat(hrs,nAnimals,1)',y,'b-','LineWidth',2,'Marker',shapes{a},'MarkerSize',markerSize)
    xticks(hrs)
    ylim([0 1])
    ylabel('DPI')
    xlabel('hours or stimulation')
    
    subplot(1,2,2);hold on
    y = meanPSS(a,hrs)';
    plot(repmat(hrs,nAnimals,1)',y,'r-','LineWidth',2,'Marker',shapes{a},'MarkerSize',markerSize)
    xticks(hrs)
    ylim([0 1])
    xlabel('hours or stimulation')
    sgtitle('mean DPI (per animal)')
end

%% Statistics

alpha = 0.05;
alphaBFC = alpha/(nHr-1);
varNames = {'hrA','hrB','p','h','bonferroni h'};
for a = 1:length(animals)
    anovaStat_V1(a) = anova1(dpiV1{a}(:,1),dpiV1{a}(:,2),'off');
    anovaStat_PSS(a) = anova1(dpiPSS{a}(:,1),dpiPSS{a}(:,2),'off');
    count = 0;

    for b = 2:nHr

        count = count+1;
        hrA(count,1) = 1;
        hrB(count,1) = b;

        x = dpiV1{a}(dpiV1{a}(:,2)==hrA(count),1);
        y = dpiV1{a}(dpiV1{a}(:,2)==hrB(count),1);
        if isempty(x) || isempty(y)
            rankP_V1(count,a) = nan;
        else
            rankP_V1(count,a) = ranksum(x,y);
        end


        x = dpiPSS{a}(dpiPSS{a}(:,2)==hrA(count),1);
        y = dpiPSS{a}(dpiPSS{a}(:,2)==hrB(count),1);
        if isempty(x) || isempty(y)
            rankP_PSS(count,a) = nan;
        else
            rankP_PSS(count,a) = ranksum(x,y);
        end


    end
    stats.dpi.v1{a} = table(hrA,hrB,rankP_V1(:,a),rankP_V1(:,a)<alpha,rankP_V1(:,a)<alphaBFC,'VariableNames',varNames);
    stats.dpi.pss{a} = table(hrA,hrB,rankP_PSS(:,a),rankP_PSS(:,a)<alpha,rankP_PSS(:,a)<alphaBFC,'VariableNames',varNames);
end
