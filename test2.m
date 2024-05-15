clear all
close all

if ispc
%     dataFold = 'D:\data'; 
%     dataFold = 'C:\Users\brand\Documents\data';
%     dataFold = 'F:\Brandon\data';
    dataFold = 'F:\Brandon\VSS2024\data';
elseif ismac
%     dataFold = '/Volumes/Lab drive/Brandon/data';
    dataFold = '/Volumes/Lab drive/Brandon/VSS2024/data';
%     dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end
physDir = fullfile(dataFold,'Ephys');
dsDir = fullfile(dataFold,'dataSets');





plr = 0;
figure; hold on
for age = 1:3

    if age == 1
        load(fullfile(dsDir,'P33to34MergeSUDat.mat'))
    elseif age == 2
        load(fullfile(dsDir,'P36to37MergeSUDat.mat'))
    elseif age == 3
        load(fullfile(dsDir,'adultMergeSUDat.mat'))
    end
    idx = dat.gu(1,:)|dat.gu(2,:);
    dat.x = dat.x(:,idx);
    dat.y = dat.y(:,idx);

    for cnd = 1:2
        

        if cnd == 1
            clr = 'k';
        elseif cnd == 2
            clr = 'c';
        end

        for u = 1:size(dat.x,2)
            x = dat.x{cnd,u};
            y = dat.y{cnd,u};
            y(y<0) = 0;
            if isnan(x) | isnan(y) | length(y)<16 
                continue
            end
            [X{age,cnd}(u,:),Y{age,cnd}(u,:),~] = alignDirTuning(x,y);
        end

        subplot(2,3,age); hold on
        plot(mean(X{age,cnd}),mean(Y{age,cnd}),clr)
        plot(repmat(mean(X{age,cnd}),2,1) , mean(Y{age,cnd})+([1;-1]*(std(Y{age,cnd})/sqrt(size(Y{age,cnd},1)))) , clr)
        ylim([-1 10])
        xticks([-180 -90 0 90 180])

        subplot(2,3,age+3); hold on
        plot( mean(X{age,cnd}) , mean(Y{age,cnd})/max(mean(Y{age,cnd})) , clr)
        ylim([-0.2 1])
        xticks([-180 -90 0 90 180])

        
    
    end
end







figure; hold on
for age = 1:3
    subplot(1,3,age);hold on
    if age == 1
        load(fullfile(dsDir,'P33to34MergeSUDat.mat'))
        clr = [0.9290 0.6940 0.1250]; %yellow
        linstyl = ':';
    elseif age == 2
        load(fullfile(dsDir,'P36to37MergeSUDat.mat'))
        clr = [0.8500 0.3250 0.0980]; %orange
        linstyl = '--';
    elseif age == 3
        load(fullfile(dsDir,'adultMergeSUDat.mat'))
        clr = [0.6350 0.0780 0.1840]; %red
        linstyl = '-';
    end
    dist{age} = 1-dat.dcv;
    dist{age}(dist{age}>1)=1;

%     %SCATTER PLOT
%     plot(dist{age}(1,:),dist{age}(2,:),'.','Color',clr,'MarkerSize',10)
%     ub = max(dist{age},[],'all');
%     plot([0 ub],[0 ub],'k--')

    %BOX PLOT
    plot(dist{age}(1:2,:),'k')
    boxplot(dist{age}(1:2,:)')
    ylim([0 1])
    xlim([0 3])


%     %CDF
%     cdf = cdfplot(dist{age}(1,:));
%     cdf.Color = clr;
%     cdf.LineStyle = '-';
%     cdf = cdfplot(dist{age}(2,:));
%     cdf.Color = clr;
%     cdf.LineStyle = '--';

end








