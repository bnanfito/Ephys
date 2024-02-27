clear
close all

% dataFold = 'D:\OneDrive - Johns Hopkins\Documents\data';
dataFold = 'D:\data';
% dataFold = 'F:\Brandon\data\sf5';
animalList = {'FEAS9','FEAT2','FEAU5','FEAU8','FEAU9','FEAV0','FEAV1'};
count = 0;
for a = 1:length(animalList)
    animal = animalList{a};
    if length(animal)~=5
        continue
    end
    
    exptList = dir(fullfile(dataFold,'Ephys',animal));
    for e = 1:height(exptList)
        exptName = exptList(e).name;
        if length(exptName)~=14 || (~strcmp(exptName(1:2),'fe') && ~strcmp(exptName(1:2),'FE'))
            continue
        end

        load(fullfile(dataFold,'Ephys',animal,exptName,[exptName '_id.mat']))
        load(fullfile(dataFold,'Ephys',animal,exptName,[exptName '.analyzer']),'-mat')
        predelay = getparam('predelay',Analyzer);
        stimTime = getparam('stim_time',Analyzer);
        postdelay = getparam('postdelay',Analyzer);
        for p = 1:length(id.probes)
            count = count+1;
            s(count).id = [id.exptId ' p' num2str(p)];
            s(count).refrTime = id.MUextractSpikes.settings{p}.refrTime;
            s(count).refrCross = id.MUextractSpikes.settings{p}.refrCross;
            s(count).spikeSamples = id.MUextractSpikes.settings{p}.spikeSamples;
            s(count).spikeRadius = id.MUextractSpikes.settings{p}.spikeRadius;
            s(count).offsetSamples = id.MUextractSpikes.settings{p}.offsetSamples;
            s(count).legacyFlag = id.MUextractSpikes.settings{p}.legacyFlag;
            s(count).predelay = predelay;
            s(count).stimTime = stimTime;
            s(count).postdelay = postdelay;
%             if isfield(id.MUthreshold,'settings')
%                 s(count).scaleFactor = id.MUthreshold.settings{p}.scaleFactor;
%             else
                load(fullfile(dataFold,'Ephys',animal,exptName,[exptName '_p' num2str(p) '_MUthreshold.mat']))
                s(count).scaleFactor = MUthresholding.threshlevel;
%             end
            load(fullfile(dataFold,'Ephys',animal,exptName,[exptName '_p' num2str(p) '_MUspkMerge.mat']))
            s(count).nSpks = length(MUspkMerge.spktimes);

        end          
        clear id
    end

end
