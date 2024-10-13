%metAna_timcourse
clear all
close all

% dataFold = '/Volumes/Lab drive/Brandon/data';
dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
anaMode = 'MU';
load(fullfile(dataFold,'dataSets','training','Train_V1Cool',anaMode,'ranksum & rPref above 2',['Train_V1Cool_' anaMode 'dataSet.mat']))

tBinSize = 0.01;
tBins = -1:tBinSize:2;
for i = 1:4

    if i == 1
        tbl = data.v1bf;
        clr = 'b';
        linStyl = '--';
    elseif i == 2
        tbl = data.v1af;
        clr = 'b';
        linStyl = '-';
    elseif i == 3
        tbl = data.pssbf;
        clr = 'r';
        linStyl = '--';
    elseif i == 4
        tbl = data.pssaf;
        clr = 'r';
        linStyl = '-';
    end

    tbl = tbl(tbl.goodUnit,:);
    nU = height(tbl);
    for u = 1:nU
        cPrefIdx = find(tbl.cndKey{u}(:,strcmp(tbl.paramKey{u},'ori'))==tbl.oriPref(u));
        if length(cPrefIdx)==2
            sizes = tbl.cndKey{u}(cPrefIdx,contains(tbl.paramKey{u},'size'));
            cPrefIdx = cPrefIdx(sizes == min(sizes));
        end
        prefTrials = tbl(u,:).fr.trialNum(:,cPrefIdx);
        prefTrialsIdx = ismember(tbl.spkTimes{u}(2,:),prefTrials);
        curSpks = tbl.spkTimes{u}(1,prefTrialsIdx); %spike times of the current unit that occur on trials of the preferred condition
        psth{i}(u,:) =  (histcounts(curSpks,tBins)/length(prefTrials))*(1/tBinSize); %multiply by inverse bin size to convert to Hz

        for s = 1:length(curSpks)
            g{i,u}(s,:) = normpdf(tBins,curSpks(s),0.01);
        end
        sdf{i}(u,:) = sum(g{i,u})/length(prefTrials);
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
        clr = 'b';
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
        clr = 'r';
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

end
legend(P,legLbl)

for i = 1:4

    if i == 1
        subplot(4,2,5);
        title('v1 before');
    elseif i == 2
        subplot(4,2,7);
        title('v1 after');
    elseif i == 3
        subplot(4,2,6);
        title('pss before');
    elseif i == 4
        subplot(4,2,8);
        title('pss after');
    end

    imagesc(sdf{i})

end
sgtitle([anaMode ' sdf, pref trials'])
