

animal = 'febm6';
unit = {'002','002','002'};
expt = {'004','021','038'};


% animal = 'febm7';
% unit = {'000','000','000'};
% expt = {'000','019','032'};

% animal = 'febm8';
% unit = {'000','000','000'};
% expt = {'010','022','032'};


stimSize = 'hemi';
anaMode = 'MU';

for e = 1:length(expt)
    [sumStats] = anaOri(animal,unit{e},expt{e},'PSS',anaMode,'D:\data',0,0);
    [PSS(e)] = plotSum(sumStats,stimSize,0);
end

for e = 1:length(expt)
    [sumStats] = anaOri(animal,unit{e},expt{e},'V1',anaMode,'D:\data',0,0);
    [V1(e)] = plotSum(sumStats,stimSize,0);
end


figure('Position',[100 100 1500 500]); hold on;

subplot(1,3,1); hold on;
title([])
    cdf = cdfplot(PSS(1).rPref);
    cdf.LineStyle = '--';
    cdf.Color = 'r';
    cdf = cdfplot(PSS(2).rPref);
    cdf.LineStyle = ':';
    cdf.Color = 'r';
    cdf = cdfplot(PSS(3).rPref);
    cdf.LineStyle = '-';
    cdf.Color = 'r';
    cdf = cdfplot(V1(1).rPref);
    cdf.LineStyle = '--';
    cdf.Color = 'b';
    cdf = cdfplot(V1(2).rPref);
    cdf.LineStyle = ':';
    cdf.Color = 'b';
    cdf = cdfplot(V1(3).rPref);
    cdf.LineStyle = '-';
    cdf.Color = 'b';
    xlabel('rPref');
    ylabel('percentile');
    legend({'PSS bf','PSS mid','PSS af','V1 bf','V1 mid','V1 af'})


subplot(1,3,2); hold on;
title([])
    cdf = cdfplot(PSS(1).dsi);
    cdf.LineStyle = '--';
    cdf.Color = 'r';
    cdf = cdfplot(PSS(2).dsi);
    cdf.LineStyle = ':';
    cdf.Color = 'r';
    cdf = cdfplot(PSS(3).dsi);
    cdf.LineStyle = '-';
    cdf.Color = 'r';
    cdf = cdfplot(V1(1).dsi);
    cdf.LineStyle = '--';
    cdf.Color = 'b';
    cdf = cdfplot(V1(2).dsi);
    cdf.LineStyle = ':';
    cdf.Color = 'b';
    cdf = cdfplot(V1(3).dsi);
    cdf.LineStyle = '-';
    cdf.Color = 'b';
    xlabel('dsi');
    xlim([0 1])
    ylabel('percentile');

subplot(1,3,3); hold on;
title([])
    cdf = cdfplot(PSS(1).ldr);
    cdf.LineStyle = '--';
    cdf.Color = 'r';
    cdf = cdfplot(PSS(2).ldr);
    cdf.LineStyle = ':';
    cdf.Color = 'r';
    cdf = cdfplot(PSS(3).ldr);
    cdf.LineStyle = '-';
    cdf.Color = 'r';
    cdf = cdfplot(V1(1).ldr);
    cdf.LineStyle = '--';
    cdf.Color = 'b';
    cdf = cdfplot(V1(2).ldr);
    cdf.LineStyle = ':';
    cdf.Color = 'b';
    cdf = cdfplot(V1(3).ldr);
    cdf.LineStyle = '-';
    cdf.Color = 'b';
    xlabel('ldr');
    xlim([0 1])
    ylabel('percentile');

sgtitle([animal ' (' anaMode ' ' stimSize ')'])






