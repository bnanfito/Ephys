
clear all
close all

animal = 'febo5';
unit = {'000','000','000'};
expt = {'005','006','008'};
probe = 'V1';
anaMode = 'MU';
dataFold = '/Volumes/Lab drive/Brandon/data';
plt = 0;
svePlt = 0;

figure; hold on
for e = 1:length(expt)

    if strcmp(expt{e},'005')
        legLbl{e} = 'no lens';
        clr{e} = 'k';
    elseif strcmp(expt{e},'006')
        legLbl{e} = '-12 lens';
        clr{e} = 'g';
    elseif strcmp(expt{e},'007')
        legLbl{e} = '+12 lens';
        clr{e} = 'r';
    elseif strcmp(expt{e},'008')
        legLbl{e} = 'opaque cover';
        clr{e} = 'c';
    end

    sumStats = anaSF(animal,unit{e},expt{e},probe,anaMode,dataFold,plt,svePlt);
    sumStats = sumStats(sumStats.goodUnit,:);

    for u = 1:height(sumStats)

        oriInd = strcmp(sumStats.paramKey{u},'ori');
        sfInd = strcmp(sumStats.paramKey{u},'s_freq');

        opIdx = squeeze(mean(sumStats.condition{u}(oriInd,:,:),2)) == sumStats.oriPref(u);
        y(u,:) = mean( sumStats.response{u}(:,:,opIdx) ,1,'omitnan');
        y(u,:) = y(u,:)./max(y(u,:));
        x = sumStats.condition{u}(sfInd,:,opIdx);

%         x = -1:0.01:2;
%         prefCond = find(sumStats.cndKey{u}(:,oriInd) == sumStats.oriPref(u) & sumStats.cndKey{u}(:,sfInd) == sumStats.sfPref(u,opIdx));
%         trialList = sumStats.fr(u).trialNum(:,prefCond);
%         tInd = ismember(sumStats.spkTimes{u}(2,:),trialList);
%         spks{u} = sumStats.spkTimes{u}(1,tInd);
%         g = normpdf(repmat(x,length(spks{u}),1)',spks{u},repmat(0.05,1,length(spks{u})));
%         sdf(u,:) = mean(g,2,'omitnan');
%         sdf_z(u,:) = ( sdf(u,:)-mean(sdf(u,:),'omitnan') ) / std(sdf(u,:),'omitnan');

    end

    sem = std(y,'omitnan')./sqrt(size(y,1));
    p(e) = plot(x,mean(y,'omitnan'),'LineWidth',2,'Color',clr{e});
    plot(repmat(x,2,1),mean(y,'omitnan')+([1;-1].*sem),'LineWidth',2,'Color',clr{e})

%     p(e) = plot(x,mean(sdf,'omitnan'),'LineWidth',2,'Color',clr{e});

%     binWidth = 0.1;
%     histBins = -1:binWidth:2;
%     p(e) = histogram([spks{:}],histBins);
%     p(e).FaceColor = 'none';
% %     p(e).FaceAlpha = 0.25;
%     p(e).EdgeColor = clr{e};
%     p(e).LineWidth = 2;
%     p(e).BinCounts = p(e).BinCounts/(height(sumStats)*length(trialList)*binWidth);

    clear spks

end

legend(p,legLbl)
xlabel('spatial freq. (cyc./vis.deg.)')
ylabel('normalized response')



