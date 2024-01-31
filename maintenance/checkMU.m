%Check MU

clear all
close all

if ispc
    dataFold = 'D:\OneDrive - Johns Hopkins\Documents\data\Ephys';
elseif ismac
    dataFold = '/Users/brandonnanfito/Library/CloudStorage/OneDrive-JohnsHopkins/Documents/data/Ephys';
end

countHas = 0;
countNeeds = 0;
animal = 'febg7';
cd(fullfile(dataFold,animal))
expts = dir;
for e = 1:height(expts)
    exptName = expts(e).name;
    if length(exptName)<10 && ~strcmp(exptName,'f')
        continue
    end
    
    disp(exptName)
    load(fullfile(dataFold,animal,exptName,[exptName '_id.mat']))
    for p = 1:length(id.probes)
        
        if isfile(fullfile(dataFold,animal,exptName,[exptName '_p' num2str(p) '_MUspkMerge.mat']))
            countHas = countHas+1;
            hasMU{countHas,1} = [exptName '_p' num2str(p)]; 
        else
            countNeeds = countNeeds+1;
            needsMU{countNeeds,1} = [exptName '_p' num2str(p)];
        end

    end

end
