% mergeExpts
clear all; close all;
if ispc
    dataFold = 'D:\data'; 
%     dataFold = 'F:\Brandon\data';
elseif ismac
%     dataFold = '/Volumes/Lab drive/Brandon/data';
    dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end
physDir = fullfile(dataFold,'Ephys');
animal = 'febj8';
units = {'003','003','003','003','003','003','003','003','003','003','003','003','003'};
expts = {'002','003','004','005','006','008','009','010','016','017','018','019','020'};
mergeID = [];
files{1,length(expts)} = [];
for e = 1:length(expts)
    mergeID = [mergeID units{e} expts{e}];
    files{e} = ['u' units{e} '_' expts{e}];
end
mergeIntan(physDir,animal,files,100,mergeID,'BRN')