clear
close all

oneDrive = 'D:\OneDrive - Johns Hopkins\Documents\data';
dataFold = 'D:\data';
animalList = {'FEAO4','FEAN6','FEAS6','FEAT1','FEAQ5'};
for a = 1:length(animalList)
    animal = animalList{a};
    if length(animal)~=5
        continue
    end
    
    exptList = dir(fullfile(oneDrive,'SummaryStats',animal));
    for e = 1:height(exptList)
        exptName = exptList(e).name;
        if length(exptName)~=14
            continue
        end

        load(fullfile(dataFold,'Ephys',animal,exptName,[exptName '_id.mat']))
        for p = 1:length(id.probes)
            load(fullfile(dataFold,'SummaryStats',animal,exptName,[exptName '_p' num2str(p) '_sumStatsMU.mat']))
            aSumStats = sumStats;
            aSpks = spks;
            clear sumStats spks
            load(fullfile(oneDrive,'SummaryStats',animal,exptName,[exptName '_p' num2str(p) '_sumStatsMU.mat']))
            bSumStats = sumStats;
            bSpks = spks;
            clear sumStats spks
            if length(unique(aSpks(1).times == bSpks(1).times)) == 1 && unique(aSpks(1).times == bSpks(1).times) == 1
                disp([exptName ' p' num2str(p) ': oneDrive = local'])
            else
                disp([exptName ' p' num2str(p) ': oneDrive ~= local'])
            end
        end          
        
    end

end
