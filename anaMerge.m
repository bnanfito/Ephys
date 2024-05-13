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

anaMode = 'MU';
for i = 1:2
    figure; hold on
    if i == 1
        load(fullfile(dataFold,'dataSets',['juviMerge' anaMode 'Dat.mat']))
        clr = [0.9290 0.6940 0.1250]; %yellow
    elseif i == 2
        load(fullfile(dataFold,'dataSets',['adultMerge' anaMode 'Dat.mat']))
        clr = [0.8500 0.3250 0.0980]; %orange
%     elseif i == 3
%         load(fullfile(dataFold,'dataSets',['adultMerge' anaMode 'Dat.mat']))
%         clr = [0.6350 0.0780 0.1840]; %red
    end
    fileID = vertcat(dat.tbl.fileID{:});
    condID = fileID(:,2);
    idxA = dat.tbl.goodUnit == 1 & condID == 1;
    idxB = dat.tbl.goodUnit == 1 & condID == 2;

    A = 1-dat.tbl(idxA,:).DCV;
    B = 1-dat.tbl(idxB,:).DCV;

    A(A>1) = 1;
    B(B>1) = 1;

    cdfA = cdfplot(A);
    cdfA.Color = 'k';
    
    cdfB = cdfplot(B);
    cdfB.Color = 'c';

    clear dat
end