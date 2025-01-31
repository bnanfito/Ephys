clear all
close all

dataFold = 'Y:\Brandon\data';
physDir = fullfile(dataFold,'Ephys');
cd(physDir)

animals = dir;
for a = 1:height(animals)
    animalName = animals(a).name;
    if length(animalName)~=5 || ~( strcmp(animalName(1:2),'FE') || strcmp(animalName(1:2),'fe') ) || ~animals(a).isdir
        continue
    end
    cd(fullfile(physDir,animalName))

    expts = dir;
    for e = 1:height(expts)
        exptName = expts(e).name;
        if length(exptName)~=14 || ~( strcmp(exptName(1:2),'FE') || strcmp(exptName(1:2),'fe') ) || ~expts(e).isdir
            continue
        end
        cd(fullfile(physDir,animalName,exptName))
        disp(exptName)
        load(fullfile(physDir,animalName,exptName,[exptName '_id.mat']))
        for p = 1:length(id.probes)
            MUspkMergeFile = fullfile(physDir,animalName,exptName,[exptName '_p' num2str(p) '_MUspkMerge.mat']);
            MUthreshFile = fullfile(physDir,animalName,exptName,[exptName '_p' num2str(p) '_MUthreshold.mat']);
            if isfile(MUspkMergeFile)
                load(MUspkMergeFile)
                load(MUthreshFile)
                if ~isfield(MUspkMerge,'info')
                    disp('updating MUspkMerge')
                    MUspkMerge.info.expname = exptName;
                    MUspkMerge.info.probeID = p;
                    MUspkMerge.info.startID = 0;
                    MUspkMerge.info.stopID = length(id.MUextractSpikes.jobEdges{p})-2;
                    MUspkMerge.info.date = id.MUextractSpikes.date{p};
                    MUspkMerge.info.jobEdges = id.MUextractSpikes.jobEdges{p};
                    MUspkMerge.info.threshold.date = id.MUthreshold.date{p};
                    MUspkMerge.info.threshold.name = id.MUthreshold.name{p};
                    MUspkMerge.info.threshold.threshlevel = MUthresholding.threshlevel;
                    MUspkMerge.info.threshold.threshlength = MUthresholding.threshlength;
                    save(MUspkMergeFile,'MUspkMerge')
                end
            end
            clear MUspkMerge MUthresholding MUspkMergeFile MUthreshFile
        end
        clear id
    end
end

disp('done')

