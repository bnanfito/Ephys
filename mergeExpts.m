% mergeExpts
clear all; close all;
if ispc
%     dataFold = 'C:\Users\brand\Documents\data';
%     dataFold = 'D:\data'; 
    dataFold = 'F:\Brandon\data';
elseif ismac
    dataFold = '/Volumes/Lab drive/Brandon/data';
%     dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end
physDir = fullfile(dataFold,'Ephys');
animal = 'febl0';
units = {'001','001','001'};
expts = {'006','010','016'};
mergeID = [];
files{1,length(expts)} = [];
for e = 1:length(expts)
    mergeID = [mergeID units{e} expts{e}];
    files{e} = ['u' units{e} '_' expts{e}];
end
mergeIntan(physDir,animal,files,100,mergeID,'BRN')