clear all
close all

anaMode = 'MU';
proj = {'V1cool_ori'};

for p = 1:length(proj)
dataFold = 'Y:\Brandon\data';
% dataFold = ['Y:\Brandon\data\dataSets\training\' proj{p} '\' anaMode '\threshold4'];
projectTbl=getProjectFiles(proj{p},1,'age','recSite','priorMFlag','priorDescr','duringMFlag','manipDescr','manipDetail');

for e = 1:height(projectTbl)

    animal = projectTbl.experimentId{e};
    unit = projectTbl.unitNr{e};
    expt = projectTbl.experimentNr{e};
    probe = projectTbl.probeId(e);
    exptName = [animal '_u' unit '_' expt];
    % necessary files
    sourceDir = fullfile(dataFold,'Ephys',animal,exptName);
    trialInfoFile = [exptName '_trialInfo.mat'];
    idFile = [exptName '_id.mat'];
    analyzerFile = [exptName '.analyzer'];
    if strcmp(anaMode,'MU')
        MUextSpkFile = [exptName '_p' num2str(probe) '_MUextractSpk.mat'];
        MUextSpkPropFile = [exptName '_p' num2str(probe) '_MUextractSpkProp.mat'];
        MUspkMergeFile = [exptName '_p' num2str(probe) '_MUspkMerge.mat'];
        MUthreshFile = [exptName '_p' num2str(probe) '_MUthreshold.mat'];
    end

%% preprocess
%     preprocessData(animal,unit,expt,probe,anaMode,0,dataFold)

%% copy to destination
%     destDir = 'Y:\Brandon\data\threshold5';
%     aniDir = fullfile(destDir,'Ephys',animal);
%     if ~isfolder(aniDir)
%         mkdir(aniDir)
%     end
%     exptDir = fullfile(aniDir,exptName);
%     if ~isfolder(exptDir)
%         mkdir(exptDir)
%     end
%     cd(exptDir)
%     if ~isfile(trialInfoFile)
%         [status(e,1),msg{e,1}] = copyfile(fullfile(sourceDir,trialInfoFile));
%     end
%     if ~isfile(idFile)
%         [status(e,2),msg{e,2}] = copyfile(fullfile(sourceDir,idFile));
%     end
%     if ~isfile(analyzerFile)
%         [status(e,3),msg{e,3}] = copyfile(fullfile(sourceDir,analyzerFile));
%     end
%     if ~isfile(MUextSpkFile)
%         [status(e,4),msg{e,4}] = copyfile(fullfile(sourceDir,MUextSpkFile));
%     end
%     if ~isfile(MUextSpkPropFile)
%         [status(e,5),msg{e,5}] = copyfile(fullfile(sourceDir,MUextSpkPropFile));
%     end
%     if ~isfile(MUspkMergeFile)
%         [status(e,6),msg{e,6}] = copyfile(fullfile(sourceDir,MUspkMergeFile));
%     end
%     if ~isfile(MUthreshFile)
%         [status(e,7),msg{e,7}] = copyfile(fullfile(sourceDir,MUthreshFile));
%     end

%% check settings
    setVars = {'exptName','probeId','nJobs','refrTime','refrCross','spikeRad','spkWindowT','spkTol','spkWindowI','offsetSamples','threshLevel','date','matches?'};
    eName{e,1} = exptName;
    probeId(e,1) = probe;
    load(fullfile(sourceDir,idFile))
    if strcmp(anaMode,'MU')
        load(fullfile(sourceDir,MUspkMergeFile))
        date{e,1} = MUspkMerge.info.date;
        if isfile(fullfile(sourceDir,MUextSpkFile))
            load(fullfile(sourceDir,MUextSpkFile))
            load(fullfile(sourceDir,MUextSpkPropFile))
            load(fullfile(sourceDir,MUthreshFile))
            nJobs(e,1) = length(extractSpk.jobStart);
            refrTime(e,1) = extractSpk.settings.refrTime;
            refrCross(e,1) = extractSpk.settings.refrCross;
            spikeRad(e,1) = extractSpk.settings.spikeRadius;
            spkWindowT{e,1} = extractSpkProp.settings.spkWindowT;
            spkTol(e,1) = extractSpkProp.settings.spkTol;
            spkWindowI(e,1) = extractSpkProp.settings.spkWindowI;
            offsetSamples(e,1) = MUthresholding.offsetSamples;
            threshLevel(e,1) = MUthresholding.threshlevel;
            dateMatch(e,1) = strcmp(date{e,1},extractSpk.date) & strcmp(date{e,1},extractSpkProp.date) & strcmp(date{e,1},MUthresholding.date);
        else
            nJobs(e,1) = length(id.MUextractSpikes.jobStart{probe});
            refrTime(e,1) = id.MUextractSpikes.settings{probe}.refrTime;
            refrCross(e,1) = id.MUextractSpikes.settings{probe}.refrCross;
            spikeRad(e,1) = id.MUextractSpikes.settings{probe}.spikeRadius;
            spkWindowT{e,1} = nan;
            spkTol(e,1) = nan;
            spkWindowI(e,1) = nan;
            offsetSamples(e,1) = id.MUextractSpikes.settings{probe}.offsetSamples;
            threshLevel(e,1) = id.MUthreshold.settings{probe}.scaleFactor;
            dateMatch(e,1) = strcmp(date{e,1},id.MUextractSpikes.date{probe}) & strcmp(date{e,1},id.MUextractSpikeProps.date{probe}) & strcmp(date{e,1},id.MUthreshold.date{probe});
        end
    end
    setTbl = table(eName,probeId,nJobs,refrTime,refrCross,spikeRad,spkWindowT,spkTol,spkWindowI,offsetSamples,threshLevel,date,dateMatch);

end

end




