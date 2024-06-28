%checkMerge

clear all; close all;
if ispc
%     dataFold = 'C:\Users\brand\Documents\data';
    dataFold = 'D:\data'; 
%     dataFold = 'F:\Brandon\data';
elseif ismac
    dataFold = '/Volumes/Lab drive/Brandon/data';
%     dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end
physDir = fullfile(dataFold,'Ephys');

animals = dir(physDir);

mergeExpts = {};
for a = 1:height(animals)
    animalID = animals(a).name;
    if length(animalID)~=5
        continue
    end

    expts = dir(fullfile(physDir,animalID));
    mergeIdx = contains(convertCharsToStrings({expts.name}),'MMM');
    mergeExpts = [mergeExpts;expts(mergeIdx).name];
    

end
