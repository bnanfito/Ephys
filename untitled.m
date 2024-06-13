clear all
close all

load('D:\data\dataSets\V1Cool_Ori16_biff_ranksum.mat','Data','animals','exptList','age')

metric = 'rBlank';

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
            xlm = [-5 15];
        end

        n(a,cc) = length(dist);

        cdf = cdfplot(dist);
        cdf.Color = clr;
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