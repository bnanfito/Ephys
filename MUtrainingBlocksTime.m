clear all
close all
load('/Volumes/NielsenHome2/Brandon/data/dataSets/training/Train_V1Cool/MU/Train_V1Cool_stimBlock_MUdataSet.mat')
%characterize activity in the training blocks - time course of response
%average across all channels and blocks


%add training type
trainType=repmat({'CoolV1'},height(projectTbl),1);
projectTbl=addvars(projectTbl,trainType);


%fix block numbers - need to deal with fact that one experiment has missing
%blocks
anim=unique(projectTbl.experimentId);
blockNr=ones(height(projectTbl),1);

for a=1:length(anim)
    idx1=find(strcmp(projectTbl.experimentId,anim{a}) & strcmp(projectTbl.recSite,'V1'));
    idx2=find(strcmp(projectTbl.experimentId,anim{a}) & strcmp(projectTbl.recSite,'PSS'));

    blockNr(idx1) = [1:length(idx1)];
%     if length(idx1)==16
%         blockNr(idx1)=[1:16]; %this addresses the issue with 0/1 as first unit
%     else
%         %in the missing case, blocks start from 0, so can use the
%         %experimentNr as is
%         for i=1:length(idx1)
%             blockNr(idx1(i))=str2double(projList.experimentNr{idx1(i)})+1;
%         end
%     end

    blockNr(idx2) = [1:length(idx2)];
%     if length(idx2)==16
%         blockNr(idx2)=[1:16];
%     else
%         for i=1:length(idx2)
%             blockNr(idx2(i))=str2double(projList.experimentNr{idx2(i)})+1;
%         end
%     end
end
projectTbl=addvars(projectTbl,blockNr);


%% compute response time courses for bidir experiments
%averaging across all animals and channels
%need to account for missing blocks

anim_CoolV1=unique(projectTbl.experimentId(strcmp(projectTbl.trainType,'CoolV1')));

binWidth=10; %ms

binVec=[-2000:10:5000];
binVec=binVec(1:end-1);
baseIdx=binVec<0;

coolV1_TcV1=zeros(1,length(binVec));
coolV1_TcPSS=zeros(1,length(binVec));



for a=1:length(anim_CoolV1)

    if strcmp(anim_CoolV1{a},'febl7') || strcmp(anim_CoolV1{a},'febm7')
        blockId=[1:14];
    else
        blockId=[1:16];
    end

    countV1=0;
    countPSS=0;

    tmpV1=zeros(1,length(binVec));
    tmpPSS=zeros(1,length(binVec));

    for b=1:length(blockId)

        pV1=find(strcmp(projectTbl.experimentId,anim_CoolV1{a}) & projectTbl.blockNr==blockId(b) & strcmp(projectTbl.recSite,'V1'));
        pPSS=find(strcmp(projectTbl.experimentId,anim_CoolV1{a}) & projectTbl.blockNr==blockId(b) & strcmp(projectTbl.recSite,'PSS'));

        exptName = projectTbl.fileBase{pV1};
        if strcmp(exptName,'febk7_u000_009') || strcmp(exptName,'febl7_u000_031')
            continue
        end

        %V1
        if ~isempty(pV1)
            MUThreshTrialData('/Volumes/NielsenHome2/Brandon/data/Ephys',exptName(1:5),exptName(8:10),exptName(12:14),projectTbl.probeId(pV1),'id',3,2,5,0)
            MUFile=fullfile('/Volumes/NielsenHome2/Brandon/data/Ephys',projectTbl.experimentId{pV1},projectTbl.fileBase{pV1},...
                [projectTbl.filePhys{pV1} '_MUThreshTrial.mat']);
            load(MUFile); %generates structures MUThresh and MUThreshInfo

            for c=1:length(MUThresh)
                [psth,~]=computePSTH(MUThresh(c).spktimes,2,5,binWidth);
                psth=psth-mean(psth(baseIdx)); %baseline correct
                tmpV1=tmpV1+psth/length(MUThresh);
                countV1=countV1+1;
            end
            clear MUThresh MUThreshInfo;
        end

        %PSS
        if ~isempty(pPSS)
            MUThreshTrialData('/Volumes/NielsenHome2/Brandon/data/Ephys',exptName(1:5),exptName(8:10),exptName(12:14),projectTbl.probeId(pPSS),'id',3,2,5,0)
            MUFile=fullfile('/Volumes/NielsenHome2/Brandon/data/Ephys',projectTbl.experimentId{pPSS},projectTbl.fileBase{pPSS},...
                [projectTbl.filePhys{pPSS} '_MUThreshTrial.mat']);
            load(MUFile); %generates structures MUThresh and MUThreshInfo

            for c=1:length(MUThresh)
                [psth,~]=computePSTH(MUThresh(c).spktimes,2,5,binWidth);
                psth=psth-mean(psth(baseIdx)); %baseline correct
                tmpPSS=tmpPSS+psth/length(MUThresh);
                countPSS=countPSS+1;
            end
            clear MUThresh MUThreshInfo;
        end
    end
    coolV1_TcV1=coolV1_TcV1+tmpV1/countV1;
    coolV1_TcPSS=coolV1_TcPSS+tmpPSS/countPSS;
    
end
coolV1_TcV1=coolV1_TcV1/length(anim_CoolV1);
coolV1_TcPSS=coolV1_TcPSS/length(anim_CoolV1);



%% plot

figure; hold on
plot(binVec,coolV1_TcV1,'b-')
plot(binVec,coolV1_TcPSS,'r-')

