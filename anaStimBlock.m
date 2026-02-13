clear all
close all
% dataFold = '/Volumes/NielsenHome2/Brandon/data';
dataFold = 'Y:\Brandon\data';
load([dataFold,'/dataSets/training/Train_V1Cool/MU/Train_V1Cool_stimBlock_MUdataSet.mat'])
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


area = 'PSS';
anaMode = 'MU';
binSize = 0.01;
bins = -1:binSize:6;

blocks = unique(projectTbl.blockNr);
figure; hold on
for b = blocks'

    idx = projectTbl.blockNr == b & strcmp(projectTbl.recSite,area);

    d{b} = vertcat(projectTbl(idx,:).sumStats{:});
    d{b} = d{b}(screenUnits(d{b},anaMode),:);
    n(b) = height(d{b});

    for i = 1:height(d{b})
    
        nTrials = numel(d{b}.fr(i).trialNum);
%         sdf = compute_sdf(d{b}(i,:));
%         sdfY{b}(i,:) = sdf.mean;
%         sdfErr{b}(i,:) = sdf.sem;
        spks = d{b}.spkTimes{i}(1,:);
        psth{b}(i,:) = histcounts(spks,bins)/(nTrials*binSize);
        psth{b}(i,:) = psth{b}(i,:)-mean(psth{b}(i,bins(2:end)<0));
    
    end
    meanPsth(b,:) = mean(psth{b});
    plot3(bins(2:end),repmat(b,1,length(bins)-1),meanPsth(b,:))
    lbl{b} = ['block #' num2str(b) '; n=' num2str(n(b))];
end
legend(lbl)
view(3)

figure;
imagesc(meanPsth)
colorbar
set(gca,'YDir','reverse')
tics = xticks*binSize;
for x = 1:length(tics)
xlbl{x} = num2str(tics(x)-1);
end
xticklabels(xlbl) 

figure;
plot(bins(2:end),mean(meanPsth))

figure
subplot(2,2,1);hold on%figure; hold on
clrs = [(1/8:1/8:1)' zeros(8,1) zeros(8,1)];
for hr = 1:8
    blks = [2*hr-1,2*hr];
    psth_hr{hr} = vertcat(psth{blks});
    meanPsth_hr(hr,:) = mean(psth_hr{hr});
%     plot3(bins(2:end),repmat(hr,1,length(bins)-1),meanPsth_hr(hr,:))
    plot(bins(2:end),meanPsth_hr(hr,:)+((hr-1)*10),'Color',clrs(hr,:))
    lbl_hr{hr} = ['hour #' num2str(hr) '; n=' num2str(size(psth_hr{hr},1))];
end
fill3([0 5 5 0],[0 0 8 8],[0 0 0 0],'k','EdgeColor','k','FaceAlpha',0)
legend(lbl_hr)
xlabel('time (s)')
ylabel('hour')
zlabel('firing rate (Hz)')
% view(3)

subplot(2,2,2);hold on;%figure;
imagesc(meanPsth_hr)
colorbar
set(gca,'YDir','reverse')
xticks([-1:6]/binSize)
tics = xticks*binSize;
for x = 1:length(tics)
xlbl{x} = num2str(tics(x)-1);
end
xticklabels(xlbl)
xlabel('time (s)')
ylabel('hour')
axis tight

binsTrain = bins;
meanTrainPsth = meanPsth;
clear projectTbl proj blocks b d n psth meanPsth lbl
load([dataFold,'/dataSets/training/Train_V1Cool/MU/threshold5/ranksum/Train_V1Cool_MUdataSet.mat'])


area = 'PSS';
anaMode = 'MU';
binSize = 0.01;
bins = -1:binSize:2;

subplot(2,2,3);hold on;%figure; hold on
for b = 1:2
    if b == 1 && strcmp(area,'V1')
        d = data.v1bf;
        linStyl = 'k-';
        lbl{b} = 'V1 before';
    elseif b == 1 && strcmp(area,'PSS')
        d = data.pssbf;
        linStyl = 'k-';
        lbl{b} = 'PSS before';
    elseif b == 2 && strcmp(area,'V1')
        d = data.v1af;
        linStyl = 'b-';
        lbl{b} = 'V1 after';
    elseif b == 2 && strcmp(area,'PSS')
        d = data.pssaf;
        linStyl = 'r-';
        lbl{b} = 'PSS after';
    end
    d = d(screenUnits(d,anaMode),:);

    for i = 1:height(d)

        nTrials = numel(d.fr(i).trialNum);
        spks = d.spkTimes{i}(1,:);
%         prefCond = find(d.condition{i}(strcmp(d.paramKey{i},'ori'),:)==d.oriPref(i));
%         prefTrials = d.fr(i).trialNum(:,prefCond);
%         spks = d.spkTimes{i}(1, ismember(d.spkTimes{i}(2,:),prefTrials) );
        psth{b}(i,:) = histcounts(spks,bins)/(nTrials*binSize);
        psth{b}(i,:) = psth{b}(i,:)-mean(psth{b}(i,bins(2:end)<0));

    end
    meanPsth(b,:) = mean(psth{b});
    plot(bins(2:end),meanPsth(b,:),linStyl)
    lbl{b} = [lbl{b} '; n=' num2str(height(d))];
end
plot(binsTrain(2:end),mean(meanTrainPsth),'-','Color',[0.5 0 0])
patch([0 1 1 0],[0 0 10 10],'k','EdgeColor','none','FaceAlpha',0.2)
legend(lbl)
xlabel('time (s)')
ylabel('firing rate (Hz)')

subplot(2,2,4);%figure;
imagesc(meanPsth)
colorbar
yticks([1,2])
yticklabels({'before','after'})
xticks([-1:2]/binSize)
tics = xticks*binSize;
for x = 1:length(tics)
xlbl{x} = num2str(tics(x)-1);
end
xticklabels(xlbl)
xlabel('time (s)')
axis tight
