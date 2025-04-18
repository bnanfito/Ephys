%metAna_V1cool_cont
clear all
close all

% dataFold = '/Volumes/Lab drive/Brandon/data';
dataFold = 'F:\Brandon\data';
load(fullfile(dataFold,'dataSets','cooling','V1cool_MU_contrast','V1cool_MU_contrast_projectTbl.mat'))

pssIdx = strcmp(projectTbl.recSite,'PSS');
v1Idx = strcmp(projectTbl.recSite,'V1');
cntrlIdx = projectTbl.priorMFlag == 0 & projectTbl.duringMFlag == 0;
coolIdx = projectTbl.duringMFlag == 1;

cntrlData = projectTbl(cntrlIdx & pssIdx,:);
cntrlAnimals = unique(cntrlData.experimentId);
for a = 1:length(cntrlAnimals)
    cntrlAge(a) = unique( cntrlData( strcmp(cntrlData.experimentId,cntrlAnimals{a}) ,: ).age );
end

coolData = projectTbl(coolIdx & pssIdx,:);
coolAnimals = unique(coolData.experimentId);
for a = 1:length(coolAnimals)
    coolAge(a) = unique( coolData( strcmp(coolData.experimentId,coolAnimals{a}) ,: ).age );
end

figure;hold on
bins = min([min(coolAge) min(cntrlAge)]):max([max(coolAge) max(cntrlAge)]);
histogram(coolAge,bins)
histogram(cntrlAge,bins)


for f = 1:height(coolData)

    [coolSumStats{f,1}] = anaCon(coolData.experimentId{f},coolData.unitNr{f},coolData.experimentNr{f},'PSS','MU',dataFold,0,0);

end
coolData.sumStats = coolSumStats;


for f = 1:height(cntrlData)

    [cntrlSumStats{f,1}] = anaCon(cntrlData.experimentId{f},cntrlData.unitNr{f},cntrlData.experimentNr{f},'PSS','MU',dataFold,0,0);

end
cntrlData.sumStats = cntrlSumStats;

animals = unique([coolAnimals;cntrlAnimals]);
for a = 1:length(animals)
    curAni = animals{a};
    ages(a) = unique([coolAge(strcmp(coolAnimals,curAni)) cntrlAge(strcmp(cntrlAnimals,curAni))]);

    curAniIdx = contains(coolData.experimentId,curAni);
    if sum(curAniIdx)==0
        sumStats.cool{a} = [];
    else
        sumStats.cool{a} = vertcat(coolData.sumStats{curAniIdx});
    end

    curAniIdx = contains(cntrlData.experimentId,curAni);
    if sum(curAniIdx)==0
        sumStats.cntrl{a} = [];
    else
        sumStats.cntrl{a} = vertcat(cntrlData.sumStats{curAniIdx});
    end

end


figure;hold on
for a = 1:length(animals)
    curAni = animals{a};

    if ~isempty(sumStats.cool{a})
        coolCF = sumStats.cool{a}.cF(sumStats.cool{a}.goodUnit);
        plot(coolAge(strcmp(coolAnimals,curAni)),mean(coolCF,'omitnan'),'c.')
    end

    if ~isempty(sumStats.cntrl{a})
        cntrlCF = sumStats.cntrl{a}.cF(sumStats.cntrl{a}.goodUnit);
        plot(cntrlAge(strcmp(cntrlAnimals,curAni)),mean(cntrlCF,'omitnan'),'k.')
    end


end


for a = 1:length(animals)

    coolTbl = sumStats.cool{a};
    cntrlTbl = sumStats.cntrl{a};
    if isempty(coolTbl) || isempty(cntrlTbl)
        continue
    end
    coolTbl = coolTbl(coolTbl.goodUnit,:);
    cntrlTbl = cntrlTbl(cntrlTbl.goodUnit,:);
    
    figure;
    subplot(1,2,1);hold on
    cntrlSpkTimes = [cntrlTbl.spkTimes{:}];
    coolSpkTimes = [coolTbl.spkTimes{:}];
    bins = [-1:0.1:2];
    h1 = histogram(coolSpkTimes(1,:),bins);
    h1.FaceColor = 'c';h1.EdgeColor = 'none';
    h2 = histogram(cntrlSpkTimes(1,:),bins);
    h2.FaceColor = 'k';h2.EdgeColor = 'none';

    for u = 1:height(coolTbl)
        meanCurve.cool{a}(u,:) = mean(coolTbl.response{u}(:,:,coolTbl.oriPref(u)),'omitnan');
        cnd.cool{a}(u,:) = coolTbl.condition{u}(strcmp(coolTbl.paramKey{u},'contrast'),:,coolTbl.oriPref(u));
    end

    for u = 1:height(cntrlTbl)
        meanCurve.cntrl{a}(u,:) = mean(cntrlTbl.response{u}(:,:,cntrlTbl.oriPref(u)),'omitnan');
        cnd.cntrl{a}(u,:) = cntrlTbl.condition{u}(strcmp(cntrlTbl.paramKey{u},'contrast'),:,cntrlTbl.oriPref(u));
    end
    
    subplot(2,2,2);hold on
    plot(cnd.cool{a}',meanCurve.cool{a}','c--','LineWidth',0.5)
    plot(cnd.cntrl{a}',meanCurve.cntrl{a}','k--','LineWidth',0.5)

    subplot(2,2,4);hold on
    x = mean(cnd.cool{a});
    y = meanCurve.cool{a};
    sem = std(y)/sqrt(size(y,1));
    plot(x,mean(y),'c-o','LineWidth',2)
    plot(repmat(x,2,1),mean(y)+([-1;1]*sem),'c','LineWidth',2)
    clear x y sem

    x = mean(cnd.cntrl{a});
    y = meanCurve.cntrl{a};
    sem = std(y)/sqrt(size(y,1));
    plot(x,mean(y),'k-o','LineWidth',2)
    plot(repmat(x,2,1),mean(y)+([-1;1]*sem),'k','LineWidth',2)
    clear x y sem

    title([animals{a} '; P' num2str(ages(a))])
end





