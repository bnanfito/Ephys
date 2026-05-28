%anaV1cool_ori_SU

clear
close all

anaMode = 'SU';
proj = 'V1cool_ori';
area = 'PSS';
dataFold = 'Y:\Brandon\data';
% dataFold = '/Volumes/NielsenHome2/Brandon/data';
% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
ageGroups = {[28 32],[33 40],[41 80],[81 120]};

%% Generate Project Table

% projTbl = getProjectFiles(proj,1,'age','recSite','penNr','priorMFlag','priorDescr',...
%                                        'duringMFlag','manipDescr','manipDetail',...
%                                        'looperNameCond1','looperNameCond2',...
%                                        'looperNameCond3','looperNameCond4',...
%                                        'looperNameCond5');
% projTbl = projTbl(strcmp(projTbl.recSite,area),:);
% animals = unique(projTbl.experimentId);
% for a = 1:length(animals)
%     aniIdx = strcmp(projTbl.experimentId,animals{a});
%     ages(a) = unique(projTbl.age(aniIdx));
% end
% load(fullfile(dataFold,'dataSets','cooling','V1cool_ori','MU','V1cool_ori_MUprojTbl.mat'))
% for e = 1:height(projTbl)
%     animal = projTbl.experimentId{e};
%     unit = projTbl.unitNr{e};
%     expt = projTbl.experimentNr{e};
%     probe = projTbl.probeId(e);
%     exptName = projTbl.fileBase{e};
%     disp(['generating sumStats for ' exptName])
%     [sumStats{e,1}] = anaOri(animal,unit,expt,probe,anaMode,dataFold,0,0);
%     nGU(e,1) = sum(screenUnits(sumStats{e,1},anaMode));
% end
% projTbl.sumStats = sumStats;
% projTbl.nGU = nGU;

load(fullfile(dataFold,'dataSets','cooling',proj,anaMode,[proj '_' anaMode 'dataSet.mat']))
% load(['/Volumes/NielsenHome2/Brandon/data/dataSets/cooling/' proj '/V1cool_' anaMode '_ori_projectTbl.mat'])
% load(['/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/cooling/' proj '/' anaMode '/V1cool_ori_' anaMode 'dataSet.mat'])

%% Organize Data

coolIdx = projTbl.duringMFlag == 1 & strcmp(projTbl.manipDescr,'Cooling');
cntrlIdx = projTbl.duringMFlag == 0 & projTbl.priorMFlag == 0;
for ag = 1:length(ageGroups)
    curAGidx = find(ages>=ageGroups{ag}(1) & ages<=ageGroups{ag}(2));
    agIdx(curAGidx) = ag;
end

agesJit = double(ages);
for u = unique(agesJit)
    if sum(agesJit==u)>1
        agesJit(agesJit==u) = agesJit(agesJit==u)+(-0.2:0.4/(sum(agesJit==u)-1):0.2);
    end
end

for a = 1:length(animals)
    aniIdx = strcmp(projTbl.experimentId,animals{a});

    %if there are repeated expts in a penetration, only take the one with highest yield of good units
    cntrlPens = unique(projTbl.penNr(cntrlIdx & aniIdx));
    for p = 1:length(cntrlPens)
        penIdx = projTbl.penNr == cntrlPens(p);
        idx = find(cntrlIdx & aniIdx & penIdx);
        if length(idx) > 1 
            idx = idx( projTbl.nGU(idx) == max(projTbl.nGU(idx)) );
        end
        datCntrl = projTbl.sumStats{idx};
        goodId = screenUnits(datCntrl,anaMode);
%         datCntrl = datCntrl(goodId,:);
        dat.cntrl{a,cntrlPens(p)} = datCntrl;
        goodIdCntrl{a,cntrlPens(p)} = goodId;

    end

    coolPens = unique(projTbl.penNr(coolIdx & aniIdx));
    for p = 1:length(coolPens)
        penIdx = projTbl.penNr == coolPens(p);
        idx = find(coolIdx & aniIdx & penIdx);
        if length(idx) > 1
            idx = idx( projTbl.nGU(idx) == max(projTbl.nGU(idx)) );
        end
        datCool = projTbl.sumStats{idx};
        goodId = screenUnits(datCool,anaMode);
%         datCool = datCool(goodId,:);
        dat.cool{a,coolPens(p)} = datCool;
        goodIdCool{a,coolPens(p)} = goodId;

    end

end

%% Calculate Animal Metrics

for a = 1:length(animals)
    if ~(isempty(dat.cntrl{a,1}) || isempty(dat.cool{a,1}))

        % For metrics of response magnitude, include all MU that pass
        % criteria before cooling
        i = goodIdCntrl{a,1}; %index MU to include
        c = 0.1; %add a constant offset to log transformed response mag.

        r1 = dat.cntrl{a,1}.rPref(i);
        r1L = log10(r1+c);
        rPref.cntrl.dist{a} = r1;
        rPref.cntrl.ave(a) = mean(r1,'omitnan');
        rPref.cntrl.sem(a) = std(r1,'omitnan')/sqrt(length(r1));
        rPref.cntrl.n(a) = length(r1);

        d = dat.cntrl{a,1}.ldr(i);
        ldr.cntrl.dist{a} = d;
        ldr.cntrl.ave(a) = mean(d,'omitnan');
        ldr.cntrl.sem(a) = sem(d);
        ldr.cntrl.n(a) = length(d);

        d = dat.cntrl{a,1}.lor(i);
        lor.cntrl.dist{a} = d;
        lor.cntrl.ave(a) = mean(d,'omitnan');
        lor.cntrl.sem(a) = sem(d);
        lor.cntrl.n(a) = length(d);

        l1 = dat.cntrl{a,1}.latency(i);
        late.cntrl.dist{a} = l1;
        late.cntrl.ave(a) = mean(l1,'omitnan');
        late.cntrl.sem(a) = std(l1,'omitnan')/sqrt(length(l1));
        late.cntrl.n(a) = length(l1);

        % Also get tuning metrics for cooled condition
        i = goodIdCool{a,1}; %index MU to include

        r2 = dat.cool{a,1}.rPref(i);
        r2L = log10(r2+c);
        rPref.cool.dist{a} = r2;
        rPref.cool.ave(a) = mean(r2,'omitnan');
        rPref.cool.sem(a) = std(r2,'omitnan')/sqrt(length(r2));
        rPref.cool.n(a) = length(r2);

%         si.dist{a} = (r2-r1)./(r2+r1);
%         % si.dist{a} = (r2L-r1L)./(r2L+r1L);
%         si.ave(a) = mean(si.dist{a},'omitnan');
%         si.sem(a) = std(si.dist{a},'omitnan')/sqrt(length(si.dist{a}));
%         si.n(a) = length(si.dist{a});
% 
%         if ~(isempty(r1)||isempty(r2))
%             rPref.stats.sRank.p(a) = signrank(r1-r2);
%             si.stats.sRank.p(a) = signrank(si.dist{a});
%         else
%             rPref.stats.sRank.p(a) = nan;
%             si.stats.sRank.p(a) = nan;
%         end
        
        d = dat.cool{a,1}.ldr(i);
        ldr.cool.dist{a} = d;
        ldr.cool.ave(a) = mean(d,'omitnan');
        ldr.cool.sem(a) = sem(d);
        ldr.cool.n(a) = length(d);

        d = dat.cool{a,1}.lor(i);
        lor.cool.dist{a} = d;
        lor.cool.ave(a) = mean(d,'omitnan');
        lor.cool.sem(a) = sem(d);
        lor.cool.n(a) = length(d);

        l2 = dat.cool{a,1}.latency(i);
        late.cool.dist{a} = l2;
        late.cool.ave(a) = mean(l2,'omitnan');
        late.cool.sem(a) = std(l2,'omitnan')/sqrt(length(l2));
        late.cool.n(a) = length(l2);

%         % For metrics of response latency, the MU included must pass
%         % criteria before and during cooling, and exclude nan values
%         i = (goodIdCntrl{a,1}&goodIdCool{a,1}) & ... %index MU to include
%             ~(isnan(dat.cntrl{a,1}.latency)|isnan(dat.cool{a,1}.latency));
% 
%         dLate.dist{a} = l2-l1;
%         dLate.ave(a) = mean(dLate.dist{a},'omitnan');
%         dLate.sem(a) = std(dLate.dist{a},'omitnan')/sqrt(length(dLate.dist{a}));
%         dLate.n(a) = length(dLate.dist{a});
% 
%         if ~(isempty(l1)||isempty(l2))
%             late.stats.sRank.p(a) = signrank(l1-l2);
%             dLate.stats.sRank.p(a) = signrank(dLate.dist{a});
%         else
%             late.stats.sRank.p(a) = nan;
%             dLate.stats.sRank.p(a) = nan;
%         end

    else
        rPref.cntrl.dist{a} = [];
        rPref.cntrl.ave(a) = nan;
        rPref.cntrl.sem(a) = nan;
        rPref.cntrl.n(a) = nan;

        rPref.cool.dist{a} = [];
        rPref.cool.ave(a) = nan;
        rPref.cool.sem(a) = nan;
        rPref.cool.n(a) = nan;

%         si.dist{a} = [];
%         si.ave(a) = nan;
%         si.sem(a) = nan;
%         si.n(a) = nan;
% 
%         rPref.stats.sRank.p(a) = nan;
%         si.stats.sRank.p(a) = nan;

        late.cntrl.dist{a} = [];
        late.cntrl.ave(a) = nan;
        late.cntrl.sem(a) = nan;
        late.cntrl.n(a) = nan;

        late.cool.dist{a} = [];
        late.cool.ave(a) = nan;
        late.cool.sem(a) = nan;
        late.cool.n(a) = nan;

%         dLate.dist{a} = [];
%         dLate.ave(a) = nan;
%         dLate.sem(a) = nan;
%         dLate.n(a) = nan;
% 
%         late.stats.sRank.p(a) = nan;
%         dLate.stats.sRank.p(a) = nan;

    end
end

%% Make data table

animalData1 = table();
animalData1.id = repmat(animals,2,1);
animalData1.age = repmat(ages',2,1);
animalData1.ageJit = repmat(agesJit',2,1);
animalData1.manip = [zeros(length(animals),1);ones(length(animals),1)];
animalData1.R_mean = [rPref.cntrl.ave';rPref.cool.ave'];
animalData1.R_sem = [rPref.cntrl.sem';rPref.cool.sem'];
animalData1.R_n = [rPref.cntrl.n';rPref.cool.n'];
animalData1.L_mean = [late.cntrl.ave';late.cool.ave'];
animalData1.L_sem = [late.cntrl.sem';late.cool.sem'];
animalData1.L_n = [late.cntrl.n';late.cool.n'];

animalData2 = table();
animalData2.id = animals;
animalData2.age = ages';
animalData2.ageJit = agesJit';
% animalData2.SI_mean = si.ave';
% animalData2.SI_sem = si.sem';
% animalData2.SI_n = si.n';
% animalData2.dL_mean = dLate.ave';
% animalData2.dL_sem = dLate.sem';
% animalData2.dL_n = dLate.n';

animalData = table();
animalData.id = animals;
animalData.age = ages';
animalData.ageJit = agesJit';
animalData.Rcntrl_mean = rPref.cntrl.ave';
animalData.Rcntrl_sem = rPref.cntrl.sem';
animalData.Rcntrl_n = rPref.cntrl.n';
animalData.Rcntrl_dist = rPref.cntrl.dist';
animalData.Rcool_mean = rPref.cool.ave';
animalData.Rcool_sem = rPref.cool.sem';
animalData.Rcool_n = rPref.cool.n';
animalData.Rcool_dist = rPref.cool.dist';
% animalData.SI_mean = si.ave';
% animalData.SI_sem = si.sem';
% animalData.SI_n = si.n';
% animalData.SI_dist = si.dist';
animalData.Lcntrl_mean = late.cntrl.ave';
animalData.Lcntrl_sem = late.cntrl.sem';
animalData.Lcntrl_n = late.cntrl.n';
animalData.Lcntrl_dist = late.cntrl.dist';
animalData.Lcool_mean = late.cool.ave';
animalData.Lcool_sem = late.cool.sem';
animalData.Lcool_n = late.cool.n';
animalData.Lcool_dist = late.cool.dist';
% animalData.dL_mean = dLate.ave';
% animalData.dL_sem = dLate.sem';
% animalData.dL_n = dLate.n';
% animalData.dL_dist = dLate.dist';

% [~,sIdx] = sort(animalData.ageJit);
% animalData = animalData(sIdx,:);


unitDataCntrl = table();
cntrl = vertcat(dat.cntrl{~ismember(1:length(animals),[3,11,13]),1});
for u = 1:height(cntrl)
    id{u,1} = [cntrl.exptName{u} '_p' num2str(cntrl.probe(u)) '_' cntrl.uInfo{u} num2str(cntrl.uID(u))];
    uAge(u,1) = ages(strcmp(cntrl.exptName{u}(1:5),animals));
    uAG(u,1) = agIdx(strcmp(cntrl.exptName{u}(1:5),animals));
    manip(u,1) = 0;
    switch [num2str(uAG(u,1)) num2str(manip(u,1))]
        case '10'
            g(u,1) = 1;
        case '11'
            g(u,1) = 2;
        case '20'
            g(u,1) = 3;
        case '21'
            g(u,1) = 4;
        case '30'
            g(u,1) = 5;
        case '31'
            g(u,1) = 6;
        case '40'
            g(u,1) = 7;
        case '41'
            g(u,1) = 8;
    end
end
unitDataCntrl.id = id;
unitDataCntrl.age = uAge;
unitDataCntrl.AG = uAG;
unitDataCntrl.manip = manip;
unitDataCntrl.g = g;
unitDataCntrl.rPref = cntrl.rPref;
unitDataCntrl.latency = cntrl.latency;
unitDataCntrl.dsi = cntrl.dsi;
unitDataCntrl.ldr = cntrl.ldr;
unitDataCntrl.osi = cntrl.osi;
unitDataCntrl.lor = cntrl.lor;
unitDataCntrl.good = screenUnits(cntrl,'MU');
clear id uAge uAG manip g

unitDataCool = table();
cool = vertcat(dat.cool{~ismember(1:length(animals),[3,11,13]),1});
for u = 1:height(cool)
    id{u,1} = [cool.exptName{u} '_p' num2str(cool.probe(u)) '_' cool.uInfo{u} num2str(cool.uID(u))];
    uAge(u,1) = ages(strcmp(cool.exptName{u}(1:5),animals));
    uAG(u,1) = agIdx(strcmp(cool.exptName{u}(1:5),animals));
    manip(u,1) = 1;
    switch [num2str(uAG(u,1)) num2str(manip(u,1))]
        case '10'
            g(u,1) = 1;
        case '11'
            g(u,1) = 2;
        case '20'
            g(u,1) = 3;
        case '21'
            g(u,1) = 4;
        case '30'
            g(u,1) = 5;
        case '31'
            g(u,1) = 6;
        case '40'
            g(u,1) = 7;
        case '41'
            g(u,1) = 8;
    end
end
unitDataCool.id = id;
unitDataCool.age = uAge;
unitDataCool.AG = uAG;
unitDataCool.manip = manip;
unitDataCool.g = g;
unitDataCool.rPref = cool.rPref;
unitDataCool.latency = cool.latency;
unitDataCool.dsi = cool.dsi;
unitDataCool.ldr = cool.ldr;
unitDataCool.osi = cool.osi;
unitDataCool.lor = cool.lor;
unitDataCool.good = screenUnits(cool,'MU');
clear id uAge uAG manip g

unitData = vertcat(unitDataCntrl,unitDataCool);
unitData = unitData(unitData.good,:);
