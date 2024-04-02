clear all;close all

if ispc
    dataFold = 'D:\data\Ephys';
elseif ismac
    dataFold = '/Users/brandonnanfito/Library/CloudStorage/OneDrive-JohnsHopkins/Documents/data/Ephys';
end
countAllExp = 0;
countHas = 0;
countNeeds = 0;
dualRecOnly = 1;

animal = 'febh3';
cd(fullfile(dataFold,animal))
exptList = dir;
exptList  = {exptList.name}';

for e = 1:length(exptList)

    exptName = exptList{e};
    if unique(ismember('fe',exptName)) | unique(ismember('FE',exptName))
        
        if length(exptName)~=14
            continue
        elseif ~isfile(fullfile(dataFold,animal,exptName,[exptName '_id.mat']))
            disp([exptName ' has no id file'])
            continue
        end

        disp(exptName)
        cd(fullfile(dataFold,animal,exptName))
        load(fullfile(dataFold,animal,exptName,[exptName '_id.mat']))

        if ~isfield(id,'probes')
            disp([exptName ' id file has no probes field'])
            continue
        elseif length(id.probes)<2 && dualRecOnly == 1 
            disp([exptName ' id file only has 1 probe'])
            continue
        end

        countAllExp = countAllExp+1;
        fileID{countAllExp,1} = exptName;
        
        nProbes = length(id.probes);
        for p = 1:nProbes
            
            fileExist(countAllExp,p) = isfile(fullfile(dataFold,animal,exptName,[exptName '_p' num2str(p) '_spkSort.mat']));
            if fileExist(countAllExp,p)
                countHas = countHas+1;
                hasSpkSort{countHas,1} = [exptName '_p' num2str(p) '_spkSort.mat'];
            else
                countNeeds = countNeeds+1;
                needsSpkSort{countNeeds,1} = [exptName '_p' num2str(p) '_spkSort.mat'];
            end
    
        end
        
    end

end

bothSorted = fileID(fileExist(:,1)&fileExist(:,2));