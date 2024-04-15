% mergeExpts
clear all; close all;
dataFold = 'D:\data';
physDir = fullfile(dataFold,'Ephys');
units = {'000','000','000','000','000','000'};
expts = {'002','006','007','010','011','014'};
mergeID = [];
files{1,length(expts)} = [];
for e = 1:length(expts)
    mergeID = [mergeID units{e} expts{e}];
    files{e} = ['u' units{e} '_' expts{e}];
end
mergeIntan(physDir,'febg9',files,100,mergeID,'BRN')