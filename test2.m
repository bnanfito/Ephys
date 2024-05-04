clear all
close all

if ispc
    dataFold = 'D:\data'; 
%     dataFold = 'C:\Users\brand\Documents\data';
%     dataFold = 'F:\Brandon\data';
elseif ismac
    dataFold = '/Volumes/Lab drive/Brandon/data';
%     dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end
physDir = fullfile(dataFold,'Ephys');

figure; hold on
for age = 1:3
    if age == 1
        load(fullfile(physDir,'P33to34MergeSUDat.mat'))
        clr = [0.9290 0.6940 0.1250]; %yellow
    elseif age == 2
        load(fullfile(physDir,'P36to37MergeSUDat.mat'))
        clr = [0.8500 0.3250 0.0980]; %orange
    elseif age == 3
        load(fullfile(physDir,'adultMergeSUDat.mat'))
        clr = [0.6350 0.0780 0.1840]; %red
    end
    dist{age} = dat.dsi;
    dist{age}(dist{age}>1)=1;
    plot(dist{age}(1,:),dist{age}(2,:),'.','Color',clr,'MarkerSize',10)

end