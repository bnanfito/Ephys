%Written by Brandon Nanfito

%This script is meant to locally process files from a project set of
%experiments by copying them from z to a local folder, then copy the
%resulting matlab files into a folder defined in 'projPath'

clear all
close all

proj = 'DSdev';
anaMode = 'MU';

dataFold = 'Y:\Brandon\data';
projPath = fullfile(dataFold,'dataSets',proj);
zPath = 'Z:\EphysNew';
localPath = 'C:\Users\brand\Documents\data';
projTbl=getProjectFiles(proj,1,'age','recSite','priorMFlag','priorDescr','duringMFlag','manipDescr','manipDetail');

for e = 1:height(projTbl)
    curAni = projTbl.experimentId{e};
    curExp = projTbl.fileBase{e};
    curProbe = projTbl.probeId(e);
    inFold = fullfile(localPath,'Ephys',curAni,curExp);
    outFold = fullfile(projPath,'Ephys',curAni,curExp);

    if ~isfolder(inFold)
        mkdir(inFold)
    end
    cd(inFold)
    copyfile(fullfile(zPath,'data',curAni,curExp))
    copyfile(fullfile(zPath,'analyzer',curAni,[curExp '.analyzer']))
    copyfile(fullfile(zPath,'processedSpikes',curAni,curExp,[curExp '_id.mat']))

    preprocessData(curAni,projTbl.unitNr{e},projTbl.experimentNr{e},curProbe,'MU',0,localPath)

    if ~isfolder(outFold)
        mkdir(outFold)
    end
    cd(outFold)
    copyfile(fullfile(inFold,[curExp '_p' num2str(curProbe) '_MUspkMerge.mat']))

    cd(fullfile(localPath,'Ephys'))
    rmdir(fullfile(localPath,'Ephys',curAni),'s')

end