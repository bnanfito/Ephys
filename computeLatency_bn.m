%computeLatency_bn   

function [late] = computeLatency_bn(data,dataInfo,trialExclude)


%time vector
binWidth=10; % bin width in ms
startBin=ceil(-dataInfo.baseTime*1000/binWidth)*binWidth; %need multiple of binWidth to make 0 an edge
stopBin=floor(dataInfo.stimTime*1000/binWidth)*binWidth;
binVec=startBin:binWidth:stopBin;


%if necessary, determine average response for different
%conditions (to determine best for every channel)
if ~isfield(data,'avgRateCond')
    %get lists of conditions
    condId=unique(dataInfo.triallist);
    %remove blank from options
    condId=condId(condId~=dataInfo.blankId);
    for i=1:length(condId)
        %find trials
        idx=find(dataInfo.triallist==condId(i) & trialExclude~=1);

        for u=1:length(data)
            data(u).avgRateCond(i)=mean(data(u).stimFrate(idx)-data(u).baseFrate(idx));
        end
    end
end

%loop through units/channels and determine latency
late=NaN(length(data),1);
cumResp=zeros(length(data),length(binVec)-1);
stdBase=NaN(length(data),1);
for u=1:length(data)
    %figure out which trials
%     switch app.ConditionsSwitch.Value
%         case 'All'
%             trialidx=find(dataInfo.triallist~=dataInfo.blankId & trialExclude~=1);
%         case 'Preferred'
            [~,maxCond]=max(data(u).avgRateCond); %channel specific preferred condition
            trialidx=find(dataInfo.triallist==maxCond & trialExclude~=1);
%     end
    NTrials=length(trialidx);

    %compute timecourse for best condition
    spkTimes{u}={};
    N=zeros(NTrials,length(binVec)-1);
    for i=1:NTrials
        spkTimes{u}{i}=data(u).spktimes{trialidx(i)};
        N(i,:)=histcounts(spkTimes{u}{i},binVec);
    end

    %mean and std values for bins for baseline
    avgN=mean(N,1)/(binWidth/1000); %in spikes/s
    avgBase=mean(avgN(binVec<0));
    stdBase(u)=std(avgN(binVec<0));

    %cum sum and its derivative
    cumN=cumsum(avgN-avgBase); %cumsum(1): value 1 in input; length(cumsum)=length(input)
    cumResp(u,:)=cumN; %copy for later
    diffN=diff(cumN); %diff(1): difference between values 2 and 1 in the input; length = length(input)-1
    cumN=cumN(1:end-1); %reduce by one to match diff

    %latency: cumsum over threshold, followed by N increasing
    %bins; example for 2 bins: cumN>thresh + diffN>0 + diffN(2:end)>0 
    %(no offset needed for first diffN since diffN(1) contains
    %difference between cumN(1) and cumN(2); next one needs to
    %be shifted to get one more increasing bin
    %to keep all matrices the same, we use circshift
    %to generate all of the shifted copies
    %also add that time has to be larger than 0 to clean up
    NShift=2;
    diffNshift=zeros(NShift,length(diffN));
    for i=0:NShift-1 %first one is not shifted, see above
        diffNshift(i+1,:)=circshift(diffN,-i);
    end
    latBin=find(cumN>2*stdBase(u) & all(diffNshift>0,1) & binVec(1:end-2)>0,1); %return only first instance
    
    if ~isempty(latBin)
        late(u)=binVec(latBin);
    end
    
end


% switch anaMode
%     case 'SU'
%         fname=fullfile(dataFold,'');
%         
%         SUlatency=app.respLat;
%         SUlatency(app.unitExclude==1)=NaN;
%         SUlatencyInfo.Condition=app.ConditionsSwitch.Value;
%         SUlatencyInfo.NrBins=app.NrBinSpin.Value;
%         SUlatencyInfo.BinWidth=app.BinWidthEdit.Value;
%        
%         save(fname,'SUlatency','SUlatencyInfo');
%     case 'MU'
%         fname=replace(app.FilenameEdit.Value,'ThreshTrial','Latency');
%         
%         MUlatency=app.respLat;
%         MUlatency(app.unitExclude==1)=NaN;
%         MUlatencyInfo.Condition=app.ConditionsSwitch.Value;
%         MUlatencyInfo.NrBins=app.NrBinSpin.Value;
%         MUlatencyInfo.BinWidth=app.BinWidthEdit.Value;
%         MUlatencyInfo.ExcludeOutliers=app.ExclOutlierSel.Value;
% 
%         save(fname,'MUlatency','MUlatencyInfo');                   
end



