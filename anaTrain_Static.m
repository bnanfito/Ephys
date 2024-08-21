%anaTrain_BiDir

clear all
close all

% load('/Volumes/Lab drive/Brandon/data/dataSets/training/Train_BiDir/projectTbl.mat')

projectTbl=getProjectFiles('Train_Static',1,'age','recSite','priorMFlag','priorDescr','duringMFlag','manipDescr','manipDetail');
dataFold = 'Y:\Brandon\data';

for e = 1:height(projectTbl)
    animal = projectTbl.experimentId{e};
    unit = projectTbl.unitNr{e};
    expt = projectTbl.experimentNr{e};
    probe = projectTbl.probeId(e);
    exptName = projectTbl.fileBase{e};
    exptDir = fullfile(dataFold,'Ephys',animal,exptName);
    disp(['generating sumStats for ' exptName])
    [sumStats{e,1}] = anaOri(animal,unit,expt,probe,'MU',dataFold,0,0);
end
projectTbl.sumStats = sumStats;