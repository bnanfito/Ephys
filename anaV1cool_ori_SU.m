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
        i = goodIdCntrl{a,1}; %index SU to include
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
        i = goodIdCool{a,1}; %index SU to include

        r2 = dat.cool{a,1}.rPref(i);
        r2L = log10(r2+c);
        rPref.cool.dist{a} = r2;
        rPref.cool.ave(a) = mean(r2,'omitnan');
        rPref.cool.sem(a) = std(r2,'omitnan')/sqrt(length(r2));
        rPref.cool.n(a) = length(r2);
        
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

    else
        rPref.cntrl.dist{a} = [];
        rPref.cntrl.ave(a) = nan;
        rPref.cntrl.sem(a) = nan;
        rPref.cntrl.n(a) = nan;

        late.cntrl.dist{a} = [];
        late.cntrl.ave(a) = nan;
        late.cntrl.sem(a) = nan;
        late.cntrl.n(a) = nan;

        ldr.cntrl.dist{a} = [];
        ldr.cntrl.ave(a) = nan;
        ldr.cntrl.sem(a) = nan;
        ldr.cntrl.n(a) = nan;

        lor.cntrl.dist{a} = [];
        lor.cntrl.ave(a) = nan;
        lor.cntrl.sem(a) = nan;
        lor.cntrl.n(a) = nan;

        rPref.cool.dist{a} = [];
        rPref.cool.ave(a) = nan;
        rPref.cool.sem(a) = nan;
        rPref.cool.n(a) = nan;

        late.cool.dist{a} = [];
        late.cool.ave(a) = nan;
        late.cool.sem(a) = nan;
        late.cool.n(a) = nan;

        ldr.cool.dist{a} = [];
        ldr.cool.ave(a) = nan;
        ldr.cool.sem(a) = nan;
        ldr.cool.n(a) = nan;

        lor.cool.dist{a} = [];
        lor.cool.ave(a) = nan;
        lor.cool.sem(a) = nan;
        lor.cool.n(a) = nan;

    end
end

%% Make data table

animalData = table();
animalData.id = animals;
% animalData.age = ages';
animalData.ageJit = agesJit';

animalData.Rcntrl_mean = rPref.cntrl.ave';
animalData.Rcntrl_sem = rPref.cntrl.sem';
animalData.Rcntrl_n = rPref.cntrl.n';
% animalData.Rcntrl_dist = rPref.cntrl.dist';

animalData.Rcool_mean = rPref.cool.ave';
animalData.Rcool_sem = rPref.cool.sem';
animalData.Rcool_n = rPref.cool.n';
% animalData.Rcool_dist = rPref.cool.dist';

animalData.LDRcntrl_mean = ldr.cntrl.ave';
animalData.LDRcntrl_sem = ldr.cntrl.sem';
animalData.LDRcntrl_n = ldr.cntrl.n';
% animalData.LDRcntrl_dist = ldr.cntrl.dist';

animalData.LDRcool_mean = ldr.cool.ave';
animalData.LDRcool_sem = ldr.cool.sem';
animalData.LDRcool_n = ldr.cool.n';
% animalData.LDRcool_dist = ldr.cool.dist';

animalData.LORcntrl_mean = lor.cntrl.ave';
animalData.LORcntrl_sem = lor.cntrl.sem';
animalData.LORcntrl_n = lor.cntrl.n';
% animalData.LORcntrl_dist = lor.cntrl.dist';

animalData.LORcool_mean = lor.cool.ave';
animalData.LORcool_sem = lor.cool.sem';
animalData.LORcool_n = lor.cool.n';
% animalData.LORcool_dist = lor.cool.dist';

animalData.Lcntrl_mean = late.cntrl.ave';
animalData.Lcntrl_sem = late.cntrl.sem';
animalData.Lcntrl_n = late.cntrl.n';
% animalData.Lcntrl_dist = late.cntrl.dist';

animalData.Lcool_mean = late.cool.ave';
animalData.Lcool_sem = late.cool.sem';
animalData.Lcool_n = late.cool.n';
% animalData.Lcool_dist = late.cool.dist';

[~,sIdx] = sort(animalData.ageJit);
animalData = animalData(sIdx,:);


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
unitDataCntrl.good = screenUnits(cntrl,anaMode);
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
unitDataCool.good = screenUnits(cool,anaMode);
clear id uAge uAG manip g

unitData = vertcat(unitDataCntrl,unitDataCool);
unitData = unitData(unitData.good,:);

%% Plot

figure;
tiledlayout(1,2)
nexttile;hold on
ageIdx = animalData.age<40;
idx = animalData.LDRcntrl_n>=5;
errorbar(animalData.ageJit(idx&ageIdx),animalData.LDRcntrl_mean(idx&ageIdx),animalData.LDRcntrl_sem(idx&ageIdx),'ko','MarkerFaceColor','k')
errorbar(animalData.ageJit(~idx&ageIdx),animalData.LDRcntrl_mean(~idx&ageIdx),animalData.LDRcntrl_sem(~idx&ageIdx),'ko')
idx = animalData.LDRcool_n>=5;
errorbar(animalData.ageJit(idx&ageIdx),animalData.LDRcool_mean(idx&ageIdx),animalData.LDRcool_sem(idx&ageIdx),'co','MarkerFaceColor','c')
errorbar(animalData.ageJit(~idx&ageIdx),animalData.LDRcool_mean(~idx&ageIdx),animalData.LDRcool_sem(~idx&ageIdx),'co')
nexttile;hold on
ageIdx = animalData.age>40;
idx = animalData.LDRcntrl_n>=5;
errorbar(animalData.ageJit(idx&ageIdx),animalData.LDRcntrl_mean(idx&ageIdx),animalData.LDRcntrl_sem(idx&ageIdx),'ko','MarkerFaceColor','k')
errorbar(animalData.ageJit(~idx&ageIdx),animalData.LDRcntrl_mean(~idx&ageIdx),animalData.LDRcntrl_sem(~idx&ageIdx),'ko')
idx = animalData.LDRcool_n>=5;
errorbar(animalData.ageJit(idx&ageIdx),animalData.LDRcool_mean(idx&ageIdx),animalData.LDRcool_sem(idx&ageIdx),'co','MarkerFaceColor','c')
errorbar(animalData.ageJit(~idx&ageIdx),animalData.LDRcool_mean(~idx&ageIdx),animalData.LDRcool_sem(~idx&ageIdx),'co')

figure;
tiledlayout(1,2)
nexttile;hold on
ageIdx = animalData.age<40;
idx = animalData.LORcntrl_n>=5;
errorbar(animalData.ageJit(idx&ageIdx),animalData.LORcntrl_mean(idx&ageIdx),animalData.LORcntrl_sem(idx&ageIdx),'ko','MarkerFaceColor','k')
errorbar(animalData.ageJit(~idx&ageIdx),animalData.LORcntrl_mean(~idx&ageIdx),animalData.LORcntrl_sem(~idx&ageIdx),'ko')
idx = animalData.LORcool_n>=5;
errorbar(animalData.ageJit(idx&ageIdx),animalData.LORcool_mean(idx&ageIdx),animalData.LORcool_sem(idx&ageIdx),'co','MarkerFaceColor','c')
errorbar(animalData.ageJit(~idx&ageIdx),animalData.LORcool_mean(~idx&ageIdx),animalData.LORcool_sem(~idx&ageIdx),'co')
nexttile;hold on
ageIdx = animalData.age>40;
idx = animalData.LORcntrl_n>=5;
errorbar(animalData.ageJit(idx&ageIdx),animalData.LORcntrl_mean(idx&ageIdx),animalData.LORcntrl_sem(idx&ageIdx),'ko','MarkerFaceColor','k')
errorbar(animalData.ageJit(~idx&ageIdx),animalData.LORcntrl_mean(~idx&ageIdx),animalData.LORcntrl_sem(~idx&ageIdx),'ko')
idx = animalData.LORcool_n>=5;
errorbar(animalData.ageJit(idx&ageIdx),animalData.LORcool_mean(idx&ageIdx),animalData.LORcool_sem(idx&ageIdx),'co','MarkerFaceColor','c')
errorbar(animalData.ageJit(~idx&ageIdx),animalData.LORcool_mean(~idx&ageIdx),animalData.LORcool_sem(~idx&ageIdx),'co')

figure;hold on
plot(unitDataCntrl.age(unitDataCntrl.good),unitDataCntrl.ldr(unitDataCntrl.good),'o')
plot(unitDataCool.age(unitDataCool.good),unitDataCool.ldr(unitDataCool.good),'o')
