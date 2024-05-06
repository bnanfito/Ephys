clear all
close all

if ispc
%     dataFold = 'D:\data'; 
%     dataFold = 'C:\Users\brand\Documents\data';
    dataFold = 'F:\Brandon\data';
elseif ismac
%     dataFold = '/Volumes/Lab drive/Brandon/data';
    dataFold = '/Volumes/Lab drive/Brandon/VSS2024/data';
%     dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end
physDir = fullfile(dataFold,'Ephys');
dsDir = fullfile(dataFold,'dataSets');




load(fullfile(dataFold,'dataSets','V1 Mono Bidir Training','coolV1Train_MUdataSet.mat'))
figure;hold on
for t = 1:2
    if t==1
        tbl = v1bf{end};
        clr = 'b';
        linStl = '--';
    elseif t==2
        tbl = v1af{end};
        clr = 'b';
        linStl = '-';
    elseif t==3
        tbl = pssbf{end};
        clr = 'r';
        linStl = '--';
    elseif t==4
        tbl = pssaf{end};
        clr = 'r';
        linStl = '-';
    end
    tbl = tbl(tbl.goodUnit==1,:);
    
    x = mean(vertcat(tbl.tuningX{:}));
    y = vertcat(tbl.tuningY{:});
    for u = 1:size(y,1)
        y(u,:) = y(u,:)./max(y(u,:));
    end
    sem = std(y)./sqrt(size(y,1));
    plot(x,mean(y),'Color',clr,'LineStyle',linStl)
    plot(repmat(x,2,1),mean(y)+([-1;1]*sem),'Color',clr)
    xticks([-180 -90 0 90 180])
    ylim([0 1])

end







% figure; hold on
% for age = 1:3
%     subplot(1,3,age);hold on
%     if age == 1
%         load(fullfile(dsDir,'P33to34MergeSUDat.mat'))
%         clr = [0.9290 0.6940 0.1250]; %yellow
%         linstyl = ':';
%     elseif age == 2
%         load(fullfile(dsDir,'P36to37MergeSUDat.mat'))
%         clr = [0.8500 0.3250 0.0980]; %orange
%         linstyl = '--';
%     elseif age == 3
%         load(fullfile(dsDir,'adultMergeSUDat.mat'))
%         clr = [0.6350 0.0780 0.1840]; %red
%         linstyl = '-';
%     end
%     dist{age} = 1-dat.dcv;
%     dist{age}(dist{age}>1)=1;
% 
% %     %SCATTER PLOT
% %     plot(dist{age}(1,:),dist{age}(2,:),'.','Color',clr,'MarkerSize',10)
% %     ub = max(dist{age},[],'all');
% %     plot([0 ub],[0 ub],'k--')
% 
%     %BOX PLOT
%     plot(dist{age}(1:2,:),'k')
%     boxplot(dist{age}(1:2,:)')
%     ylim([0 1])
%     xlim([0 3])
% 
% 
% %     %CDF
% %     cdf = cdfplot(dist{age}(1,:));
% %     cdf.Color = clr;
% %     cdf.LineStyle = '-';
% %     cdf = cdfplot(dist{age}(2,:));
% %     cdf.Color = clr;
% %     cdf.LineStyle = '--';
% 
% end











