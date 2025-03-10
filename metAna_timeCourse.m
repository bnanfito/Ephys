%metAna_timcourse
clear all
close all

% dataFold = '/Volumes/Lab drive/Brandon/data';
% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
dataFold = 'Y:\Brandon\data';
anaMode = 'MU';
proj = 'Train_V1Cool';
allTrials = 1;
load(fullfile(dataFold,'dataSets','training',proj,anaMode,[proj '_' anaMode 'dataSet.mat']))

tBinSize = 0.01;
tBins = -1:tBinSize:2;
for i = 1:4

    if i == 1
        tbl = data.v1bf;
    elseif i == 2
        tbl = data.v1af;
    elseif i == 3
        tbl = data.pssbf;
    elseif i == 4
        tbl = data.pssaf;
    end

    tbl = tbl(tbl.goodUnit,:);
    nU = height(tbl);
    for u = 1:nU

        oriIdx = strcmp(tbl.paramKey{u},'ori');
        sizeIdx = contains(tbl.paramKey{u},'size');
        cPrefIdx = find(tbl.cndKey{u}(:,oriIdx)==tbl.oriPref(u));
        if sum(sizeIdx)>0
            cPrefSizes = tbl.cndKey{u}(cPrefIdx,sizeIdx);
            cPrefIdx = cPrefIdx(cPrefSizes == min(cPrefSizes));
            hemiIdx = tbl.cndKey{u}(:,sizeIdx) == min(cPrefSizes);
        end

        if allTrials == 1
            if sum(sizeIdx)>0
                hemiTrials = tbl(u,:).fr.trialNum(:,hemiIdx);
                hemiTrials = hemiTrials(:);
                hemiSpkIdx = ismember(tbl.spkTimes{u}(2,:),hemiTrials);
                curSpks = tbl.spkTimes{u}(1,hemiSpkIdx);
                nTrials = length(unique(tbl.spkTimes{u}(2,hemiSpkIdx)));
            else
                curSpks = tbl.spkTimes{u}(1,:);
                nTrials = length(unique(tbl.spkTimes{u}(2,:)));
            end
        else
            prefTrials = tbl(u,:).fr.trialNum(:,cPrefIdx);
            prefTrials = prefTrials(:);
            prefSpkIdx = ismember(tbl.spkTimes{u}(2,:),prefTrials);
            curSpks = tbl.spkTimes{u}(1,prefSpkIdx); %spike times of the current unit that occur on trials of the preferred condition
            nTrials = length(prefTrials);
        end
        psth{i}(u,:) =  (histcounts(curSpks,tBins)/nTrials)*(1/tBinSize); %multiply by inverse bin size to convert to Hz

        for s = 1:length(curSpks)
            g{i,u}(s,:) = normpdf(tBins,curSpks(s),0.01);
        end
        sdf{i}(u,:) = sum(g{i,u})/nTrials;
        late{i}(u) = tbl.latency(u);

    end
    [~,lateIdx{i}] = sort(late{i});
    sdf{i} = sdf{i}(lateIdx{i},:);


end


%% Plotting

figure; hold on
kernel = normpdf(-3:tBinSize*50:3);
for i = 1:4

    if i == 1
        tbl = data.v1bf;
        clr = 'b';
        linStyl = '--';
        legLbl{i} = 'v1 before';
        subplot(2,2,1);hold on
    elseif i == 2
        tbl = data.v1af;
        clr = 'c';
        linStyl = '-';
        legLbl{i} = 'v1 after';
        subplot(2,2,1);hold on
    elseif i == 3
        tbl = data.pssbf;
        clr = 'r';
        linStyl = '--';
        legLbl{i} = 'pss before';
        subplot(2,2,2);hold on
    elseif i == 4
        tbl = data.pssaf;
        clr = 'm';
        linStyl = '-';
        legLbl{i} = 'pss after';
        subplot(2,2,2);hold on
    end

%     x = tBins(2:end);
%     y = mean(psth{i},'omitnan');
%     sem = std(psth{i},'omitnan')/sqrt(size(psth{i},1));
%     patch([x fliplr(x)],[conv(y+sem,kernel,'same') fliplr(conv(y-sem,kernel,'same'))],clr,'EdgeColor','none','FaceAlpha',0.2)
%     plot(x,conv(y,kernel,'same'),[clr linStyl],'LineWidth',1)

    x = tBins;
    y = mean(sdf{i},'omitnan');
    sem = std(sdf{i},'omitnan')/sqrt(size(sdf{i},1));
    patch([x fliplr(x)],[y+sem fliplr(y-sem)],clr,'EdgeColor','none','FaceAlpha',0.2)
    P(i) = plot(x,y,[clr linStyl],'LineWidth',1);
    xlabel('time (s) rel. to stim.')
    ylabel('SDF')

end
legend(P,legLbl)

for i = 1:4

    if i == 1
        subplot(4,2,5);
        ttl = 'v1 before';
    elseif i == 2
        subplot(4,2,7);
        ttl = 'v1 after';
    elseif i == 3
        subplot(4,2,6);
        ttl = 'pss before';
    elseif i == 4
        subplot(4,2,8);
        ttl = 'pss after';
    end

    imagesc(sdf{i})
    title(ttl)

end
if allTrials == 1
    sgtitle([proj ' ' anaMode ' sdf, all trials'])
else
    sgtitle([proj ' ' anaMode ' sdf, pref trials'])
end
