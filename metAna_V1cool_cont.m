%metAna_V1cool_cont
clear all
close all

load('/Volumes/Lab drive/Brandon/data/dataSets/cooling/V1cool_MU_contrast/V1cool_MU_contrast_projectTbl.mat');

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

    [coolSumStats{f,1}] = anaCon(coolData.experimentId{f},coolData.unitNr{f},coolData.experimentNr{f},'PSS','MU','/Volumes/Lab drive/Brandon/data',0,0);

end
coolData.sumStats = coolSumStats;


for f = 1:height(cntrlData)

    [cntrlSumStats{f,1}] = anaCon(cntrlData.experimentId{f},cntrlData.unitNr{f},cntrlData.experimentNr{f},'PSS','MU','/Volumes/Lab drive/Brandon/data',0,0);

end
cntrlData.sumStats = cntrlSumStats;

animals = unique([coolAnimals;cntrlAnimals]);
for a = 1:length(animals)
    curAni = animals{a};

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
        plot(coolAge(strcmp(coolAnimals,curAni)),mean(coolCF,'omitnan'),'co')
    end

    if ~isempty(sumStats.cntrl{a})
        cntrlCF = sumStats.cntrl{a}.cF(sumStats.cntrl{a}.goodUnit);
        plot(cntrlAge(strcmp(cntrlAnimals,curAni)),mean(cntrlCF,'omitnan'),'ko')
    end


end





