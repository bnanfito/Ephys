clear all
close all

if ispc
%     dataFold = 'D:\data'; 
%     dataFold = 'C:\Users\brand\Documents\data';
    dataFold = 'F:\Brandon\data';
elseif ismac
%     dataFold = '/Volumes/Lab drive/Brandon/data';
    dataFold = '/Volumes/Lab drive/Brandon/VSS2024/data';
%     dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end
physDir = fullfile(dataFold,'Ephys');
dsDir = fullfile(dataFold,'dataSets');

figure; hold on
for age = 1:3
    subplot(1,3,age);hold on
    if age == 1
        load(fullfile(dsDir,'P33to34MergeSUDat.mat'))
        clr = [0.9290 0.6940 0.1250]; %yellow
        linstyl = ':';
    elseif age == 2
        load(fullfile(dsDir,'P36to37MergeSUDat.mat'))
        clr = [0.8500 0.3250 0.0980]; %orange
        linstyl = '--';
    elseif age == 3
        load(fullfile(dsDir,'adultMergeSUDat.mat'))
        clr = [0.6350 0.0780 0.1840]; %red
        linstyl = '-';
    end
    dist{age} = 1-dat.dcv;
    dist{age}(dist{age}>1)=1;

%     %SCATTER PLOT
%     plot(dist{age}(1,:),dist{age}(2,:),'.','Color',clr,'MarkerSize',10)
%     ub = max(dist{age},[],'all');
%     plot([0 ub],[0 ub],'k--')

    %BOX PLOT
    boxplot(dist{age}(:,:)')

%     %CDF
%     cdf = cdfplot(dist{age}(1,:));
%     cdf.Color = clr;
%     cdf.LineStyle = '-';
%     cdf = cdfplot(dist{age}(2,:));
%     cdf.Color = clr;
%     cdf.LineStyle = '--';

end



