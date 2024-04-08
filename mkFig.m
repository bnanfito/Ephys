
clear
close all

anaMode = 'SU';
% animal = 'febk8';
% expts = {[animal '_u000_000'],[animal '_u002_016'],[animal '_u002_033']};
animal = 'febk7';
expts = {[animal '_u000_001'],[animal '_u000_025'],[animal '_u000_042']};

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
figure;hold on
countCDF = 0;
for e = 1:length(expts)
    load(fullfile(physDir,animal,expts{e},[expts{e} '_id.mat']))
    for p = 1:length(id.probes)
        fileName = fullfile(sumDir,animal,expts{e},[expts{e} '_p' num2str(p) '_sumStats' anaMode '.mat']);
        load(fileName)
        sumStats.dsi(sumStats.dsi>1)=1;
        sumStats = sumStats(sumStats.goodUnit==1,:);
        if strcmp(id.probes(p).area,'V1')
            v1{e} = sumStats;
            area = 'V1';
            clr = 'b';
        elseif strcmp(id.probes(p).area,'PSS')
            pss{e} = sumStats;
            area = 'PSS';
            clr = 'r';
        end
        if ~isempty(sumStats)
        cdf = cdfplot(sumStats.dsi);
        countCDF = countCDF+1;
        cdf.LineStyle = linStyl{e};
        cdf.LineWidth = 2;
        cdf.Color = clr;
        lbl{countCDF} = [area ' ' exptLbls{e} '; n=' num2str(height(sumStats))];
        else

        end

        clear sumStats spks
    end
    legend(lbl,'Location','southeast')

end
title(animal)
ylabel('percentile')
xlabel('DSI')



