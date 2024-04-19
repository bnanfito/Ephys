%compute latencies across files

projList=getProjectFiles('DSdevMU',1,'age','recSite');

%% loop through files
muData=table;

binWidth=10; %ms
%all MU files have 1sec baseline and 1sec stim period
startBin=ceil(-1000/binWidth)*binWidth; %need multiple of binWidth to make 0 an edge
stopBin=floor(1000/binWidth)*binWidth;
binVec=[startBin:binWidth:stopBin];


for p=1:height(projList)
%    disp(i)
%for p=1:1

    %load file
    MUFile=fullfile('D:\ephys\DSdev',projList.experimentId{p},projList.fileBase{p},...
        [projList.filePhys{p} '_MUThreshTrial.mat']);
    load(MUFile); %generates structures MUThresh and MUThreshInfo

    trialExclude=MUThreshFlagOutlier3(MUThresh,MUThreshInfo,0);

    %load id file to get Z
    idFile=fullfile('D:\ephys\DSdev',projList.experimentId{p},projList.fileBase{p},...
        [projList.fileBase{p} '_id.mat']);
    load(idFile);

    %find sites to include
    %inclusion criteria: pref rate>X, anova across all conditions <.05
    for c=1:length(MUThresh)
        
        %include?
        pVal=anova1(MUThresh(c).stimNspk-MUThresh(c).baseNspk,MUThreshInfo.triallist,'off');

        if pVal<.05

            %compute direction response
            nDir=length(MUThreshInfo.domval);
            nRep=length(MUThreshInfo.triallist)/length(unique(MUThreshInfo.triallist));
            FRateDir=zeros(nRep,nDir);
            Ftmp=MUThresh(c).stimFrate-MUThresh(c).baseFrate;
            Ftmp(trialExclude)=NaN;
            for i=1:nDir
                idx=find(MUThreshInfo.triallist==i);
                FRateDir(:,i)=Ftmp(idx);
            end

            dirData=getDirTuning(FRateDir,MUThreshInfo.domval,0);
        
            if dirData.prefResp>0

                %prep output
                sumData=struct;

                %find correct trials
                prefDirIdx=find(MUThreshInfo.domval==dirData.prefDir);
                mxTrials=find(MUThreshInfo.triallist==prefDirIdx);
                NTrials=length(mxTrials);
                sumData.prefResp=dirData.prefResp;
     
                %compute timecourse for best condition
                spkT={};
                baseR=[];
                N=zeros(NTrials,length(binVec)-1);
                for i=1:NTrials
                    if ~trialExclude(mxTrials(i))
                        spkT{i}=MUThresh(c).spktimes{mxTrials(i)};
                        N(i,:)=histcounts(spkT{i},binVec);
                    end
                end
           
                %bin and compute latency
                %avg first
                avgN=mean(N,1)/(binWidth/1000); %in spikes/s

                %mean and std of rate before stim onset - this needs to be per
                %bin too to be reasonable
                avgBase=mean(avgN(binVec<0));
                stdBase=std(avgN(binVec<0));

                %cum sum
                cumN=cumsum(avgN-avgBase); %cumsum(1): value 1 in input
                diffN=diff(cumN); %diff(1): difference between values 2 and 1 in the input

                %criterion 1: cumsum over threshold and 2 increasing bins
                idx=find(cumN(1:end-2)>2*stdBase & diffN(1:end-1)>0 & diffN(2:end)>0 ...
                    & binVec(1:end-3)>0,1);
                if ~isempty(idx)
                    sumData.latCh=binVec(idx);
                else
                    sumData.latCh=NaN;
                end

                %criterion 2: cumsum over threshold and 3 increasing bins
                idx=find(cumN(1:end-3)>2*stdBase & ...
                    diffN(1:end-2)>0 & diffN(2:end-1)>0 & diffN(3:end)>0 ...
                    & binVec(1:end-4)>0,1);
                if ~isempty(idx)
                    sumData.latCh2=binVec(idx);
                else
                    sumData.latCh2=NaN;
                end


                %other data
                sumData.channelId=c;
                sumData.animal=projList.experimentId(p);
                sumData.age=projList.age(p);
                sumData.recSite=projList.recSite(p);

                %depth
                chId=MUThresh(c).detCh;
                sumData.zpos=id.probes(projList.probeId(i)).z(chId);

                %transfer into output
                muData=[muData;struct2table(sumData)];

            end %if resp>0
        end %if pval
    end %for channel
end %for file


%% generate 2 more tables with separate entries for V1 and PSS
V1rows=strcmp(muData.recSite,'V1');
V1Data=muData(V1rows,:);

PSSrows=strcmp(muData.recSite,'PSS');
PSSData=muData(PSSrows,:);

%save tables
save('y:\kristina\ds paper\DSdev_MULatency.mat','muData','V1Data','PSSData')




