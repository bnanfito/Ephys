%anaDsDev

%This function contructs the project table (including processed data
%tables produced by 'anaOri.m') for the DsDev project. It also copies the
%necessary files for processing the data into 'dataFold'. All subsequent
%analysis for this project is done on the DSdev_dataSet.mat file this function
%produces.

clear all
close all
proj = 'DSdev';
anaMode = 'SU';

dataFold = fullfile('Y:\Brandon\data\dataSets',proj);
zPath = 'Z:\EphysNew';
projectTbl=getProjectFiles(proj,1,'age','recSite','priorMFlag','priorDescr','duringMFlag','manipDescr','manipDetail');

for e = 1:height(projectTbl)
    animal = projectTbl.experimentId{e};
    unit = projectTbl.unitNr{e};
    expt = projectTbl.experimentNr{e};
    probe = projectTbl.probeId(e);
    exptName = projectTbl.fileBase{e};
    
    % make directories if necessary
    aniDir = fullfile(dataFold,'Ephys',animal);
    if ~isfolder(aniDir)
        mkdir(aniDir)
    end
    exptDir = fullfile(aniDir,exptName);
    if ~isfolder(exptDir)
        mkdir(exptDir)
    end
    cd(exptDir)

    % copy files if necessary
    spkSortFile = [exptName '_p' num2str(probe) '_spkSort.mat'];
    trialInfoFile = [exptName '_trialInfo.mat'];
    idFile = [exptName '_id.mat'];
    analyzerFile = [exptName '.analyzer'];
    
    if ~isfile(spkSortFile)
        [status(e,1),msg{e,1}] = copyfile(fullfile(zPath,'processedSpikes',animal,exptName,spkSortFile));
    end
    if ~isfile(trialInfoFile)
        [status(e,2),msg{e,2}] = copyfile(fullfile(zPath,'processedSpikes',animal,exptName,trialInfoFile));
    end
    if ~isfile(idFile)
        [status(e,3),msg{e,3}] = copyfile(fullfile(zPath,'processedSpikes',animal,exptName,idFile));
    end
    if ~isfile(analyzerFile)
        [status(e,4),msg{e,4}] = copyfile(fullfile(zPath,'analyzer',animal,analyzerFile));
    end

    % generate summary statistics
    disp(['generating sumStats for ' exptName])
    [sumStats{e,1}] = anaOri(animal,unit,expt,probe,anaMode,dataFold,0,0);
    sumStats{e,1} = sumStats{e,1}(sumStats{e,1}.goodUnit,:);
    
end
projectTbl.sumStats = sumStats;

%sort rows of project table by age
[~,sortIdx] = sort(projectTbl.age);
projectTbl = projectTbl(sortIdx,:);









