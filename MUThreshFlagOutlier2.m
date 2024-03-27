function [trialInclude,outTrials,sortIdx]=MUThreshFlagOutlier2(MUThresh,MUThreshInfo,plotOpt)

%flag trials as outliers; outliers based on probe-wide metric

%shortcuts
nTrial=length(MUThreshInfo.triallist);
condId=unique(MUThreshInfo.triallist);
nCond=length(condId);
nRep=length(MUThreshInfo.triallist)/nCond;
nCh=length(MUThresh);

%sort trials based on condition, then reorganize into matrix that has every
%condition in a separate column
[~,sortIdx]=sort(MUThreshInfo.triallist);
sortTrial=reshape(sortIdx,nRep,nCond);


%compute summary metric for every trial
spkMat=[MUThresh.baseNspk]; %all sites for all trials
spkMat=reshape(spkMat,nTrial,nCh);
trialSum=sum(spkMat,2);

%now reorganize based on sort
trialSum=trialSum(sortTrial);

%find outliers
% outProbe=isoutlier(trialSum,'median');
for c = 1:size(trialSum,2)
    r = trialSum(:,c);
    scale = -1/(sqrt(2)*erfcinv(3/2));
    scaledMad = scale*median(abs(r-median(r)));
    outProbe(:,c) = r>(median(r)+(scaledMad*3));
end
trialInclude=~outProbe(:); %outliers = 0
outTrials = sortTrial(outProbe);

%plot if selected
if plotOpt==1
    chAll=zeros(nRep*nCond,nCh+1);
    for c=1:nCh
        chTmp=NaN(nRep,nCond);
        for i=1:length(condId)
            chTmp(:,i)=MUThresh(c).baseNspk(sortTrial(:,i));
        end
        chAll(:,c)=chTmp(:); %cond 1, all rep; cond 2, all rep etc
    end

    idxOut=find(trialInclude==0);

    figure
    t=tiledlayout(1,2,'TileSpacing','None');
    nexttile
    imagesc(chAll);
    colormap('jet')
    xlabel('Ch')
    ylabel('Trials (sorted by cond)')
    set(gca,'YTick',[1:nRep:nTrial])
    if ~isempty(idxOut)
        hold on
        for i=1:length(idxOut)
            plot([1 nCh],[idxOut(i) idxOut(i)],'r-')
        end
    end
    
    nexttile
    plot(trialSum(:),[1:nTrial],'-o')
    set(gca,'YLim',[0.5 nTrial+0.5],'YTickLabel',[],'YTick',[1:nRep:nTrial],'YDir','reverse','YGrid','on')
    xlabel('Sum Ch')
    if ~isempty(idxOut)
        hold on
        xL=get(gca,'XLim');
        for i=1:length(idxOut)
            plot([xL(1) xL(2)],[idxOut(i) idxOut(i)],'r-')
        end
    end

end




