clear all
% close all

if ispc
    dataFold = 'D:\data'; 
%     dataFold = 'C:\Users\brand\Documents\data';
%     dataFold = 'F:\Brandon\data';
elseif ismac
%     dataFold = '/Volumes/Lab drive/Brandon/data';
    dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end

animal = 'febg2';
unit = 'MMM';
expt = '001013001016001017001020001021';
probe = 1;
U = 1;
anaMode = 'SU';

exptName = [animal '_u' unit '_' expt];
load(fullfile(dataFold,'Ephys',animal,exptName,[anaMode 'dat.mat']),'uDat')
dat = uDat(uDat.uID == U,:);


figure; 
lw = 1.5;
timeLims = [-1 2];
subplot(2,2,3); hold on
maxTrial = 0;
for f = 1:height(dat)
    nTrials = max(max(dat.fr{f}.trial));
    yOffset = nTrials*(dat.fileID{f}(1)-1);
    x = dat.raster{f}(1,:);
    y = dat.raster{f}(2,:)+yOffset;
    if dat.fileID{f}(2)==2
        patch([-1 -0.9 -0.9 -1],[maxTrial+1 maxTrial+1 maxTrial+nTrials maxTrial+nTrials],'c','EdgeColor','none')
    end
    plot(x,y,'k.')
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
    plot(x,mean(y),[clr 'o'],'LineWidth',lw,'MarkerSize',5)
    plot(repmat(x,2,1),mean(y)+([-1;1]*sem),clr,'LineWidth',lw)
    xlim([0 360])
    xlabel('Direction of Motion (deg)')
    ylabel('Firing Rate (Hz)')

end

[wvData,wvAvg,wvError]=waveformReader('D:\data\Ephys','D:\data\Ephys',animal,unit,expt,probe,U,20,0);
subplot(2,2,4);hold on
plot(wvData(:,:)','Color',[0.8 0.8 0.8])
plot(wvAvg,'k','LineWidth',lw)
plot((wvAvg+([-1;1]*wvError))','k--','LineWidth',lw)

sgtitle([animal ' ' unit ' ' expt ' p' num2str(probe) ' u' num2str(U)])

