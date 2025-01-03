%anaDsDev

proj = 'DSdev';
anaMode = 'MU';

dataFold = fullfile('Y:\Brandon\data\dataSets',proj);
projectTbl=getProjectFiles(proj,1,'age','recSite','priorMFlag','priorDescr','duringMFlag','manipDescr','manipDetail');

for e = 1:height(projectTbl)
    animal = projectTbl.experimentId{e};
    unit = projectTbl.unitNr{e};
    expt = projectTbl.experimentNr{e};
    probe = projectTbl.probeId(e);
    exptName = projectTbl.fileBase{e};
    aniDir = fullfile(dataFold,'Ephys',animal);
    exptDir = fullfile(aniDir,exptName);
%     spkSortFile = fullfile(exptDir,[exptName '_p' num2str(probe) '_spkSort.mat']);
%     analyzerFile = fullfile(exptDir,[exptName '.analyzer']);
%     trialInfo = fullfile()
    if ~isfolder(aniDir)
        mkdir(aniDir)
    end
    if ~isfolder(exptDir)
        mkdir(exptDir)
    end
%     if ~isfile()
%         copyfile()
%     end




    disp(['generating sumStats for ' exptName])
    [sumStats{e,1}] = anaOri(animal,unit,expt,probe,anaMode,dataFold,0,0);
end
projectTbl.sumStats = sumStats;











