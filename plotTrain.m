%plotTrain
clear all
close all

% dataFold = '/Volumes/Lab drive/Brandon/data';
dataFold = 'Y:\Brandon\data';
area = 'PSS';
anaMode = 'MU';
before = {'febj4_u000_003'};
during = {};
after = {'febj4_u000_024'};

%% Organize Data

for e = 1:3

    if e == 1
        list = before;
    elseif e == 2
        list = during;
    elseif e == 3
        list = after;
    end

    for l = 1:length(list)

        dat{l,e} = [];
        exptName = list{l};
        animal = exptName(1:5);
        unit = exptName(8:10);
        expt = exptName(12:14);
        [dat{l,e}] = anaOri(animal,unit,expt,area,anaMode,dataFold,0,0);
        dat{l,e} = dat{l,e}(screenUnits(dat{l,e},anaMode),:);
    end

end

%% PLOT

metric = 'ldr';

figure; hold on
count = 0;
clrs = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
for e = 1:3

    if e == 1
        list = before;
        linStyl = ':';
    elseif e == 2
        list = during;
        linStyl = '--';
    elseif e == 3
        list = after;
        linStyl = '-';
    end

    for l = 1:length(list)

        if isempty(dat{l,e})
            continue
        end
    
        count = count+1;
        p(count) = cdfplot(dat{l,e}{:,metric});
        p(count).LineStyle = linStyl;
        p(count).LineWidth = 2;
        p(count).Color = clrs{l};
        pLbl{count} = [dat{l,e}.exptName{1} '; n=' num2str(height(dat{l,e}))];

    end

end
xlim([0 1]);
xlabel(metric)
ylabel('percentile')
legend(p,pLbl)
title(area)
clear p pLbl



figure; hold on
count = 0;
for e = 1:3

    datMerge{e} = vertcat(dat{:,e});
    count = count+1;
    if e == 1
        linStyl = ':';
        pLbl{count} = ['before; n=' num2str(height(datMerge{e}))];
    elseif e == 2
        linStyl = '--';
        pLbl{count} = ['mid; n=' num2str(height(datMerge{e}))];
    elseif e == 3
        linStyl = '-';
        pLbl{count} = ['after; n=' num2str(height(datMerge{e}))];
    end
    
    p(count) = cdfplot(datMerge{e}{:,metric});
    p(count).LineStyle = linStyl;
    p(count).LineWidth = 2;
    p(count).Color = clrs{1};

end
xlim([0 1])
xlabel(metric)
ylabel('percentile')
legend(p,pLbl)
title(area)



