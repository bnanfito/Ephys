
clear
close all

anaMode = 'MU';
metric = 'DSI';
animals = {'febk7','febk8'};
ckSum = 1;

figure;
for a = 1:length(animals)

    animal = animals{a};

    if strcmp(animal,'febk7')
        expts = {[animal '_u000_001'],...
            [animal '_u000_025'],...
            [animal '_u000_042']};
    elseif strcmp(animal,'febk8')
        expts = {[animal '_u000_000'],...
            [animal '_u002_016'],....
            [animal '_u002_033']};
    end

    if ispc
        dataFold = 'D:\data';
    %     dataFold = 'D:\OneDrive - Johns Hopkins\Documents\data';
    %     dataFold = 'F:\Brandon\data';
    elseif ismac
        dataFold = '/Volumes/Lab drive/Brandon/data';
    end
    physDir = fullfile(dataFold,'Ephys');
    figDir = fullfile(dataFold,'Figures');
    sumDir = fullfile(dataFold,'SummaryStats');
    
    linStyl = {':','--','-'};
    exptLbls = {'Before Training',' Training Midpoint','After Training'};
    
    subplot(1,length(animals)+1,a);hold on
    countCDF = 0;
    for e = 1:length(expts)
        load(fullfile(physDir,animal,expts{e},[expts{e} '_id.mat']))
        for p = 1:length(id.probes)
            fileName = fullfile(sumDir,animal,expts{e},[expts{e} '_p' num2str(p) '_sumStats' anaMode '.mat']);
            if ckSum == 1 && isfile(fileName)
                load(fileName)
            else
                [sumStats,spks] = plotUnits(animal,unit,expt,probe,anaMode,visTest,alpha,plt,saveSum,saveFigs,dataFold);
            end
            sumStats.dsi(sumStats.dsi>1)=1;
            sumStats = sumStats(sumStats.goodUnit==1,:);
            if strcmp(id.probes(p).area,'V1')
                v1{a,e} = sumStats;
                areas{p} = 'V1';
                clr = 'b';
            elseif strcmp(id.probes(p).area,'PSS')
                pss{a,e} = sumStats;
                areas{p} = 'PSS';
                clr = 'r';
            end
            if ~isempty(sumStats)
                if strcmp(metric,'DCV') || strcmp(metric,'dcv')
                    cdf = cdfplot(1-sumStats.dcv);
                elseif strcmp(metric,'DSI') || strcmp(metric,'dsi')
                    cdf = cdfplot(sumStats.dsi);
                end
                countCDF = countCDF+1;
                cdf.LineStyle = linStyl{e};
                cdf.LineWidth = 2;
                cdf.Color = clr;
                lbl{countCDF} = [areas{p} ' ' exptLbls{e} '; n=' num2str(height(sumStats))];
            end
    
            clear sumStats spks
        end
    
    end
    legend(lbl,'Location','southeast')
    title(animal)
    ylabel('percentile')
    if strcmp(metric,'DCV') || strcmp(metric,'dcv')
        xlabel('1-DCV')
    elseif strcmp(metric,'DSI') || strcmp(metric,'dsi')
        xlabel('DSI')
    end

    clear lbl
end

subplot(1,length(animals)+1,length(animals)+1); hold on
countCDF = 0;
for e = 1:3
    for a = 1:length(areas)

        if a == 1
            dat = vertcat(v1{:,e});
            clr = 'b';
            area = 'V1';
        elseif a == 2
            dat = vertcat(pss{:,e});
            clr = 'r';
            area = 'PSS';
        end

        if ~isempty(dat)
            countCDF = countCDF+1;
            cdf = cdfplot(1-dat.dcv);
            cdf.Color = clr;
            cdf.LineStyle = linStyl{e};
            cdf.LineWidth = 2;
            lbl{countCDF} = [area ' ' exptLbls{e} '; n=' num2str(height(dat)) ];
        end
        
        clear dat
    end
end
xlabel('1-DCV')
ylabel('percentile')
title('all animals')
legend(lbl,'Location','southeast')
