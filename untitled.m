clear all
close all

load('/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/V1Cool_Ori16_monoChemi_ranksum.mat','Data','animals','exptList','age')

metrics = {'rPref','dsi','ldr'};

for m = 1:length(metrics)

metric = metrics{m};

figure; hold on
for a = 1:size(Data,1)
    animal = animals{a};
    subplot(4,4,a);hold on

    for cc = 1:size(Data,2)
        if isempty(Data{a,cc})
            continue
        end
        
        if cc == 1
            clr = 'k';
        elseif cc == 2
            clr = 'c';
        elseif cc == 3
            clr = 'm';
        end
        goodInd = Data{a,cc}.goodUnit==1;
        if sum(goodInd)==0
            continue
        end

        if strcmp(metric,'dsi')
            dist = Data{a,cc}.dsi(goodInd);
            dist(dist>1) = 1;
            xlbl = 'DSI';
            xlm = [0 1];
        elseif strcmp(metric,'ldr')
            dist = Data{a,cc}.ldr(goodInd);
            dist(dist>1) = 1;
            xlbl = 'Ldir';
            xlm = [0 1];
        elseif strcmp(metric,'rPref')
            dist = Data{a,cc}.rPref(goodInd);
            dist2 = Data{a,cc}.rNull(goodInd);
            xlbl = 'BCFR (Hz)';
            xlm = [0 max(dist)];
        elseif strcmp(metric,'rBlank')
            dist = mean([Data{a,cc}.rBlank{:}],'omitnan');
            xlbl = 'BCFR (Hz)';
            xlm = [0 max(dist)];
        end

        n(a,cc) = length(dist);

        cdf = cdfplot(dist);
        cdf.Color = clr;
        cdf.LineWidth = 2;
        if strcmp(metric,'rPref')
            cdf = cdfplot(dist2);
            cdf.Color = clr;
            cdf.LineStyle = '--';
%             if cc == 1
%                 legend({'rPref','rNull'})
%             end
        end
        xlim(xlm)
        xlabel(xlbl)
        ylabel('percentile')

    end
    title([animal '; P' num2str(age{a})])

end
sgtitle(metric)

end


% figure; hold on
countU = 0;
curves = [];
for a = 1

    for cc = 1:size(Data,2)

        if cc == 1
            clr = 'k';
        elseif cc == 2
            clr = 'c';
        elseif cc == 3
            clr = 'm';
        end
        
        curData = Data{a,cc};
        if isempty(curData)
           continue
        end
        goodInd = curData.goodUnit==1;
        curData = curData(goodInd,:);
        if isempty(curData)
           continue
        end
        nU(a,cc) = size(curData,1);

        for u = 1:nU(a,cc)
            countU = countU + 1;
            
            r{countU} = curData.response{u,1};
            meanR(:,countU) = mean(r{countU},'omitnan');
            c(:,countU) = curData.condition{u,1};
            uID(:,countU) = repmat(countU,size(meanR,1),1);

            [cAligned(:,countU),meanRaligned(:,countU),i] = alignDirTuning(c(:,countU)',meanR(:,countU)');
            rAligned{countU} = r{countU}(:,i);
            uIDaligned(:,countU) = uID(i,countU);

            x = cAligned(:,countU);
            y = rAligned{countU};
            yMean = meanRaligned(:,countU);
            z = uIDaligned(:,countU);

%             plot3(x,z,yMean,clr)

            figure;hold on;
            plot(x,y',[clr '.'])
            plot(x,yMean,clr,'LineWidth',2)



        end

    end
    

    

end