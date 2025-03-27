%plotTrain
clear all
close all

dataFold = '/Volumes/Lab drive/Brandon/data';
area = 'PSS';
anaMode = 'MU';
before = {'febp7_u000_001'};
during = {'febp7_u000_002'};
after = {'febp7_u000_003'};

for e = 1:3

    dat{e} = [];
    if e == 1
        list = before;
    elseif e == 2
        list = during;
    elseif e == 3
        list = after;
    end

    for l = 1:length(list)
        exptName = list{l};
        animal = exptName(1:5);
        unit = exptName(8:10);
        expt = exptName(12:14);
        [sumStats] = anaOri(animal,unit,expt,area,anaMode,dataFold,0,0);
        sumStats = sumStats(sumStats.goodUnit==1,:);
        dat{e} = vertcat(dat{e},sumStats);
    end

end

figure; hold on
for e = 1:3

    p = cdfplot(dat{e}.osi);

end

