
clear all
close all 

if ispc
    dataFold = 'D:\data'; 
%     dataFold = 'C:\Users\brand\Documents\data';
%     dataFold = 'F:\Brandon\data';
elseif ismac
    dataFold = '/Volumes/Lab drive/Brandon/data';
%     dataFold = '/Volumes/Lab drive/Brandon/VSS2024/data';
%     dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end
physDir = fullfile(dataFold,'Ephys');


% figure; hold on
% anaMode = 'SU';
% for i = 1:3
% 
%     if i == 1
%         load(fullfile(dataFold,'dataSets',['P33to34Merge' anaMode 'Dat.mat']))
%         clr = [0.9290 0.6940 0.1250]; %yellow
%     elseif i == 2
%         load(fullfile(dataFold,'dataSets',['P36to37Merge' anaMode 'Dat.mat']))
%         clr = [0.8500 0.3250 0.0980]; %orange
%     elseif i == 3
%         load(fullfile(dataFold,'dataSets',['adultMerge' anaMode 'Dat.mat']))
%         clr = [0.6350 0.0780 0.1840]; %red
%     end
%     idx = dat.gu(1,:)|dat.gu(2,:);
%     A = dat.rn(1,idx); 
%     B = dat.rn(2,idx);
% %     A(A>1) = 1; B(B>1) = 1;
%     n(i) = length(A);
% 
%     subplot(2,2,1);hold on
% %     delta = B-A;
%     delta = (B-A)./(A+B);
%     DELTA{i,1} = delta;
%     cdf = cdfplot(delta);
%     cdf.Color = clr; cdf.LineWidth = 2;
%     xlabel('delta rNull')
%     xlim([-1 1])
%     ylabel('percentile')
%     title('')
% 
%     if i == 3
%         legend({['P33-34 ' anaMode '; n=' num2str(n(1))],['P36-37 ' anaMode '; n=' num2str(n(2))],['adult ' anaMode '; n=' num2str(n(3))]},'Location','southeast')
%     end
% 
%     subplot(2,2,2);hold on
%     ub = 16;lb = -3;
%     plot(A,B,'.','Color',clr,'MarkerSize',10)
% %     plot(log2(A),log2(B),'.','Color',clr,'MarkerSize',10)
%     plot([lb ub],[lb ub],'k--')
%     xlim([lb ub])
%     ylim([lb ub])
%     xlabel('rNull before cooling (ms)')
%     ylabel('rNull after cooling (ms)')
% 
%     
%     
%     
%     A = dat.rp(1,idx);
%     B = dat.rp(2,idx);
% 
%     subplot(2,2,3);hold on
%     delta = (B-A)./(A+B);
%     DELTA{i,2} = delta;
%     cdf = cdfplot(delta);
%     cdf.Color = clr; cdf.LineWidth = 2;
%     xlabel('delta rPref')
%     xlim([-1 1])
%     ylabel('percentile')
%     title('')
%     
%     subplot(2,2,4);hold on
%     plot(A,B,'.','Color',clr,'MarkerSize',10)
% %     plot(log2(A),log2(B),'.','Color',clr,'MarkerSize',10)
%     plot([0 30],[0 30],'k--')
%     xlabel('rPref before cooling (ms)')
%     ylabel('rPref after cooling (ms)')
% 
% end










animal = 'febl6';
unit = 'MMM';
expt = '001001001002001003001004001005001006';
probe = 1;
for U = [5 7 8]

anaMode = 'SU';

exptName = [animal '_u' unit '_' expt];
load(fullfile(physDir,animal,exptName,[anaMode 'dat.mat']),'uDat')
dat = uDat(uDat.uID == U,:);


figure; 
lw = 1.5;
timeLims = [-1 2];
subplot(2,2,3); hold on
maxTrial = 0;
for f = 1:height(dat)
    nTrials = max(max(dat.fr{f}.trial));
    yOffset = nTrials*(dat.fileID{f}(1)-1);
    if isempty(dat.raster{f})
        continue
    end
    if dat.fileID{f}(2)== 1 || dat.fileID{f}(2)== 3
        clr = 'k';
    elseif dat.fileID{f}(2)== 2
%         clr = 'c';
        clr = 'k';
        patch([timeLims fliplr(timeLims)],[maxTrial+1 maxTrial+1 maxTrial+nTrials maxTrial+nTrials],'c','EdgeColor','none','FaceAlpha',0.2)
    end
    x = dat.raster{f}(1,:);
    y = dat.raster{f}(2,:)+yOffset;
    plot(x,y,[clr '.'])
    maxTrial = maxTrial+nTrials;
end
patch([0 1 1 0],[0 0 maxTrial+1 maxTrial+1],'k','EdgeColor','none','FaceAlpha',0.1)
xlim(timeLims)
xlabel('Time (sec)')
ylim([0 maxTrial+1])
ylabel('Trial #')


fileID = vertcat(dat.fileID{:});
for c = 1:2
    if c == 1 %warm
        bit = (fileID(:,2) == 1) | (fileID(:,2) == 3);
        clr = 'k';
    elseif c == 2 %cool
        bit = fileID(:,2) == 2;
        clr = 'c';
    end

    subplot(2,2,1)
    gRaster = [dat(bit,:).raster{:}];
    if isempty(gRaster)
        continue
    end
    nTrials = 0;
    for f = find(bit')
        nTrials = nTrials + max(max(dat.fr{f}.trial));
    end
    binW = 0.1;
    hist = histogram(gRaster(1,:),[-1:binW:2],'EdgeColor','none','FaceColor',clr,'FaceAlpha',0.5);hold on
    hist.BinCounts = (hist.BinCounts/nTrials)/binW;
    ylabel('Firing Rate (Hz)')
    xlim(timeLims)
    xlabel('Time (sec)')

    subplot(2,2,2); hold on
    x = mean( vertcat(dat.tuningX{bit}) );
    y = vertcat(dat.tuningY{bit});
    sem = std(y)/sqrt(size(y,1));
    plot(x,mean(y),[clr '-o'],'LineWidth',lw,'MarkerSize',5)
    plot(repmat(x,2,1),mean(y)+([-1;1]*sem),clr,'LineWidth',lw)
    xlim([0 360])
    xticks([0 90 180 270])
    xlabel('Direction of Motion (deg)')
    ylabel('Firing Rate (Hz)')

end

[wvData,wvAvg,wvError]=waveformReader(physDir,physDir,animal,unit,expt,probe,U,20,0);
subplot(2,2,4);hold on
plot(wvData(:,:)','Color',[0.8 0.8 0.8])
plot(wvAvg,'k','LineWidth',lw)
plot((wvAvg+([-1;1]*wvError))','k--','LineWidth',lw)

sgtitle([animal ' ' unit ' ' expt ' p' num2str(probe) ' u' num2str(U)])

end


