clear all
% close all

if ispc
%     dataFold = 'D:\data'; 
%     dataFold = 'C:\Users\brand\Documents\data';
%     dataFold = 'F:\Brandon\data';
dataFold = 'F:\Brandon\VSS2024\data';
elseif ismac
%     dataFold = '/Volumes/Lab drive/Brandon/data';
    dataFold = '/Volumes/Lab drive/Brandon/VSS2024/data';
%     dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end
physDir = fullfile(dataFold,'Ephys');
dsDir = fullfile(dataFold,'dataSets');

anaMode = 'SU';
for i = 1:3
    figure; hold on
    if i == 1
        load(fullfile(dataFold,'dataSets',['P33to34Merge' anaMode 'Dat.mat']))
        clr = [0.9290 0.6940 0.1250]; %yellow
        ttl = 'P33-34';
    elseif i == 2
        load(fullfile(dataFold,'dataSets',['P36to37Merge' anaMode 'Dat.mat']))
        clr = [0.8500 0.3250 0.0980]; %orange
        ttl = 'P36-37';
    elseif i == 3
        load(fullfile(dataFold,'dataSets',['adultMerge' anaMode 'Dat.mat']))
        clr = [0.6350 0.0780 0.1840]; %red
        ttl = 'adult';
    end
    fileID = vertcat(dat.tbl.fileID{:});
    condID = fileID(:,2);
    idxA = dat.tbl.goodUnit == 1 & condID == 1;
    idxB = dat.tbl.goodUnit == 1 & condID == 2;

    A = dat.tbl(idxA,:).DSI;
    B = dat.tbl(idxB,:).DSI;

    A(A>1) = 1;
    B(B>1) = 1;

    cdfA = cdfplot(A);
    cdfA.Color = 'k';
    
    cdfB = cdfplot(B);
    cdfB.Color = 'c';

    xlim([0 1])
    title(ttl)
    clear dat
end