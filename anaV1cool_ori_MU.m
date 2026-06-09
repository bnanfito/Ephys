%anaV1cool_ori_MU
clear
close all

anaMode = 'MU';
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

% Calculate Animal Metrics
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

        r2 = dat.cool{a,1}.rPref(i);
        r2L = log10(r2+c);
        rPref.cool.dist{a} = r2;
        rPref.cool.ave(a) = mean(r2,'omitnan');
        rPref.cool.sem(a) = std(r2,'omitnan')/sqrt(length(r2));
        rPref.cool.n(a) = length(r2);

        si.dist{a} = (r2-r1)./(r2+r1);
        % si.dist{a} = (r2L-r1L)./(r2L+r1L);
        si.ave(a) = mean(si.dist{a},'omitnan');
        si.sem(a) = std(si.dist{a},'omitnan')/sqrt(length(si.dist{a}));
        si.n(a) = length(si.dist{a});

        if ~(isempty(r1)||isempty(r2))
            rPref.stats.sRank.p(a) = signrank(r1-r2);
            si.stats.sRank.p(a) = signrank(si.dist{a});
        else
            rPref.stats.sRank.p(a) = nan;
            si.stats.sRank.p(a) = nan;
        end

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

        si.dist{a} = [];
        si.ave(a) = nan;
        si.sem(a) = nan;
        si.n(a) = nan;

        rPref.stats.sRank.p(a) = nan;
        si.stats.sRank.p(a) = nan;

        ldr.cntrl.dist{a} = [];
        ldr.cntrl.ave(a) = nan;
        ldr.cntrl.sem(a) = nan;
        ldr.cntrl.n(a) = nan;

        ldr.cool.dist{a} = [];
        ldr.cool.ave(a) = nan;
        ldr.cool.sem(a) = nan;
        ldr.cool.n(a) = nan;

        lor.cntrl.dist{a} = [];
        lor.cntrl.ave(a) = nan;
        lor.cntrl.sem(a) = nan;
        lor.cntrl.n(a) = nan;

        lor.cool.dist{a} = [];
        lor.cool.ave(a) = nan;
        lor.cool.sem(a) = nan;
        lor.cool.n(a) = nan;

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

        late.stats.sRank.p(a) = nan;
%         dLate.stats.sRank.p(a) = nan;

    end
end

% Make animal data table
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

animalData.SI_mean = si.ave';
animalData.SI_sem = si.sem';
animalData.SI_n = si.n';
% animalData.SI_dist = si.dist';

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

% animalData.dL_mean = dLate.ave';
% animalData.dL_sem = dLate.sem';
% animalData.dL_n = dLate.n';
% % animalData.dL_dist = dLate.dist';

[~,sIdx] = sort(animalData.ageJit);
animalData = animalData(sIdx,:);

% Make unit data table
% this will be a table of all units in order to analyze the distribution of
% each metric across individual observations (MU)
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

unitData = vertcat(unitDataCntrl,unitDataCool); %combine control and cool data into one table
unitData = unitData(unitData.good,:); %only keep units that pass inclusion criteria

%format tables in a way that prism can group easily, do not assume rows are
%paired observations
for g = 1:8
    nU(g) = sum(unitData.g==g)';
end
uIdx = 1:max(nU);

AG1control_rPref = unitData.rPref(unitData.g==1); AG1control_rPref(end+1:max(nU)) = nan;
AG1cool_rPref = unitData.rPref(unitData.g==2); AG1cool_rPref(end+1:max(nU)) = nan;
AG2control_rPref = unitData.rPref(unitData.g==3); AG2control_rPref(end+1:max(nU)) = nan;
AG2cool_rPref = unitData.rPref(unitData.g==4); AG2cool_rPref(end+1:max(nU)) = nan;
AG3control_rPref = unitData.rPref(unitData.g==5); AG3control_rPref(end+1:max(nU)) = nan;
AG3cool_rPref = unitData.rPref(unitData.g==6); AG3cool_rPref(end+1:max(nU)) = nan;
AG4control_rPref = unitData.rPref(unitData.g==7); AG4control_rPref(end+1:max(nU)) = nan;
AG4cool_rPref = unitData.rPref(unitData.g==8); AG4cool_rPref(end+1:max(nU)) = nan;
RPREF = table(AG1control_rPref,AG1cool_rPref,AG2control_rPref,AG2cool_rPref,AG3control_rPref,AG3cool_rPref,AG4control_rPref,AG4cool_rPref);

AG1control_ldr = unitData.ldr(unitData.g==1); AG1control_ldr(end+1:max(nU)) = nan;
AG1cool_ldr = unitData.ldr(unitData.g==2); AG1cool_ldr(end+1:max(nU)) = nan;
AG2control_ldr = unitData.ldr(unitData.g==3); AG2control_ldr(end+1:max(nU)) = nan;
AG2cool_ldr = unitData.ldr(unitData.g==4); AG2cool_ldr(end+1:max(nU)) = nan;
AG3control_ldr = unitData.ldr(unitData.g==5); AG3control_ldr(end+1:max(nU)) = nan;
AG3cool_ldr = unitData.ldr(unitData.g==6); AG3cool_ldr(end+1:max(nU)) = nan;
AG4control_ldr = unitData.ldr(unitData.g==7); AG4control_ldr(end+1:max(nU)) = nan;
AG4cool_ldr = unitData.ldr(unitData.g==8); AG4cool_ldr(end+1:max(nU)) = nan;

AG1control_lor = unitData.lor(unitData.g==1);
AG1cool_lor = unitData.lor(unitData.g==2);
AG2control_lor = unitData.lor(unitData.g==3);
AG2cool_lor = unitData.lor(unitData.g==4);
AG3control_lor = unitData.lor(unitData.g==5);
AG3cool_lor = unitData.lor(unitData.g==6);
AG4control_lor = unitData.lor(unitData.g==7);
AG4cool_lor = unitData.lor(unitData.g==8);

AG1control_late = unitData.latency(unitData.g==1);
AG1cool_late = unitData.latency(unitData.g==2);
AG2control_late = unitData.latency(unitData.g==3);
AG2cool_late = unitData.latency(unitData.g==4);
AG3control_late = unitData.latency(unitData.g==5);
AG3cool_late = unitData.latency(unitData.g==6);
AG4control_late = unitData.latency(unitData.g==7);
AG4cool_late = unitData.latency(unitData.g==8);

AG1control_dsi = unitData.dsi(unitData.g==1);
AG1cool_dsi = unitData.dsi(unitData.g==2);
AG2control_dsi = unitData.dsi(unitData.g==3);
AG2cool_dsi = unitData.dsi(unitData.g==4);
AG3control_dsi = unitData.dsi(unitData.g==5);
AG3cool_dsi = unitData.dsi(unitData.g==6);
AG4control_dsi = unitData.dsi(unitData.g==7);
AG4cool_dsi = unitData.dsi(unitData.g==8);

AG1control_osi = unitData.osi(unitData.g==1);
AG1cool_osi = unitData.osi(unitData.g==2);
AG2control_osi = unitData.osi(unitData.g==3);
AG2cool_osi = unitData.osi(unitData.g==4);
AG3control_osi = unitData.osi(unitData.g==5);
AG3cool_osi = unitData.osi(unitData.g==6);
AG4control_osi = unitData.osi(unitData.g==7);
AG4cool_osi = unitData.osi(unitData.g==8);

%% V1 dependent vs independent

%for each animal, compute the proportion of MU that are V1 dependent or
%independent
for a = 1:length(animals)
    if  isempty(goodIdCntrl{a,1})||isempty(goodIdCool{a,1})
        v1_ind{a,1} = [];
        v1_dep{a,1} = [];
        badU{a,1} = [];
        disInh{a,1} = [];

        v1_depProp(a) = nan;
        v1_indProp(a) = nan;
        badU_prop(a) = nan;
        continue
    end
    v1_ind{a,1} = goodIdCntrl{a,1}&goodIdCool{a,1};
    v1_dep{a,1} = goodIdCntrl{a,1}&~goodIdCool{a,1};
    badU{a,1} = ~goodIdCntrl{a,1}&~goodIdCool{a,1};
    disInh{a,1} = ~goodIdCntrl{a,1}&goodIdCool{a,1};

    v1_depProp(a) = sum(v1_dep{a})/length(v1_dep{a});
    v1_indProp(a) = sum(v1_ind{a})/length(v1_ind{a});
    badU_prop(a) = sum(badU{a})/length(badU{a});
    disInh_prop(a) = sum(disInh{a})/length(disInh{a});

end

%sort according to age
[sAge,sIdx] = sort(ages);

%plot
figure; bar([v1_indProp(sIdx)' v1_depProp(sIdx)' badU_prop(sIdx)' disInh_prop(sIdx)'],'stacked')
xticks(1:length(sAge))
xticklabels(sAge)
legend({'V1 independent','V1 dependent','bad units','disinhibited units'})

%latency as a function of si (for MU that pass before cooling)
clear r1 r2
for a = 1:length(animals)
    if isempty(dat.cntrl{a,1})||isempty(dat.cool{a,1})
        r1{a} = [];
        r2{a} = [];
        lat1{a} = [];
        SI{a} = [];
        AGE{a} = [];
        continue
    end
    r1{a} = dat.cntrl{a,1}(v1_ind{a}|v1_dep{a},:).rPref;
    r2{a} = dat.cool{a,1}(v1_ind{a}|v1_dep{a},:).rPref;
    lat1{a} = dat.cntrl{a,1}(v1_ind{a}|v1_dep{a},:).latency;
    SI{a} = (r2{a}-r1{a})./(r2{a}+r1{a});
    AGE{a} = repmat(ages(a),length(SI{a}),1);
end

%plot
figure;
for ag = 1:4
    subplot(2,2,ag)
    scatter(vertcat(SI{agIdx==ag}),vertcat(lat1{agIdx==ag}))
    title(['ag' num2str(ag)])
    xlabel('latency')
    ylabel('SI')
end

figure;
plot3(vertcat(AGE{:}),vertcat(lat1{:}),vertcat(SI{:}),'o')
xlabel('age')
ylabel('control latency')
zlabel('SI')

figure;hold on
for ag = 1:4
    scatter(vertcat(SI{agIdx==ag}),vertcat(lat1{agIdx==ag}))
    lbl{ag} = ['ag' num2str(ag) ' n=' num2str(length(vertcat(lat1{agIdx==ag})))];
    
end
legend(lbl)
clear lbl

figure;hold on
for ag = 1:4
    cdfplot(vertcat(lat1{agIdx==ag}))
end
legend({'ag1','ag2','ag3','ag4'})

%% supp fig 1: per animal scatter - control v cool

% for ag = 1:length(ageGroups)
%     nAni(ag) = sum(ages>=ageGroups{ag}(1) & ages<=ageGroups{ag}(2));
% end
% nC = max(nAni);
nC = 8;
figure;
ag = 1;
count = 0;
[~,sortIdx] = sort(ages);
datAG.cntrl{length(ageGroups)} = [];
datAG.cool{length(ageGroups)} = [];
for ani = sortIdx

    dA =  dat.cntrl{ani,1};
    dB = dat.cool{ani,1};
    if isempty(dA)|isempty(dB)
        continue
    end
    if ages(ani)>ageGroups{ag}(2)
        ag = ag+1;
        count = 1;
    else
        count = count+1;
    end
    datAG.cntrl{ag} = [datAG.cntrl{ag};dA];
    datAG.cool{ag} = [datAG.cool{ag};dB];
    subplot(length(ageGroups),nC,(nC*(ag-1))+count); hold on
    x = dA.rPref;
    y = dB.rPref;
    m = max([x;y]);
    idA = screenUnits(dA,anaMode);
    idB = screenUnits(dB,anaMode);
    plot([0 m],[0 m],'k')
    scatter(x,y,'k')
    sPlt(1) = scatter(x(idA&~idB),y(idA&~idB),'k','MarkerFaceColor','k');
    sPlt(2) = scatter(x(idA&idB),y(idA&idB),'k','MarkerFaceColor','b');
    axis tight
    axis square
    box on
    title([animals{ani} '; P' num2str(ages(ani))])
    xlabel('control')
    ylabel('cool V1')
    if ag == length(ageGroups) && ani == sortIdx(end)
        legend(sPlt,{'V1 dependent','V1 independent'},'Location','best')
    end
end

figure
for ag = 1:length(ageGroups)
    
    dA = datAG.cntrl{ag};
    dB = datAG.cool{ag};
    idA = screenUnits(dA,anaMode);
    idB = screenUnits(dB,anaMode);
    
    subplot(1,length(ageGroups),ag);hold on
    x = dA.ldr;
    y = dB.ldr;
    m = max([x;y]);
    plot([0 m],[0 m],'k')
    scatter(x,y,'k')
    scatter(x(idA&~idB),y(idA&~idB),'k','MarkerFaceColor','k')
    scatter(x(idA&idB),y(idA&idB),'k','MarkerFaceColor','b')
    xlabel('pre-cooling Ldir')
    ylabel('v1 cooled Ldir')
    axis square
    axis tight
    box on

end 


%% Plot

figure('Position',[100 100 1000 700])

a = 4; 
sp = 1;
for u = [42 52 51]
subplot(3,3,sp); hold on
d1 = dat.cntrl{a};
d2 = dat.cool{a};
nTrials1 = max(dat.cntrl{a}.fr(u).trialNum,[],'all');
nTrials2 = max(dat.cool{a}.fr(u).trialNum,[],'all');
bw = 0.01;
bins = -1:bw:2;
h1 = histcounts(d1.spkTimes{u}(1,:),bins)/(nTrials1*bw);
s1 = compute_sdf(d1(u,:));
h2 = histcounts(d2.spkTimes{u}(1,:),bins)/(nTrials2*bw);
s2 = compute_sdf(d2(u,:));
maxH = max([s1.ci95 s2.ci95],[],'all')+1;
patch([0 1 1 0],[0 0 maxH maxH],'k','EdgeColor','none','FaceAlpha',0.2)
% plot(bins(2:end),h1,'k')
% plot(bins(2:end),h2,'c')
plot(s1.time,s1.mean,'k','LineWidth',1)
patch([s1.time fliplr(s1.time)],[s1.ci95(1,:) fliplr(s1.ci95(2,:))],'k','EdgeColor','none','FaceAlpha',0.2)
plot(s2.time,s2.mean,'c','LineWidth',1)
patch([s2.time fliplr(s2.time)],[s2.ci95(1,:) fliplr(s2.ci95(2,:))],'c','EdgeColor','none','FaceAlpha',0.2)
% plot(s1.time,s1.cMean(s1.cPref,:),'k','LineWidth',1)
% patch([s1.time fliplr(s1.time)],[s1.cMean(s1.cPref,:)+s1.cSem(s1.cPref,:) fliplr(s1.cMean(s1.cPref,:)-s1.cSem(s1.cPref,:))],'k','EdgeColor','none','FaceAlpha',0.2)
% plot(s2.time,s2.cMean(s2.cPref,:),'c','LineWidth',1)
% patch([s2.time fliplr(s2.time)],[s2.cMean(s2.cPref,:)+s2.cSem(s2.cPref,:) fliplr(s2.cMean(s2.cPref,:)-s2.cSem(s2.cPref,:))],'c','EdgeColor','none','FaceAlpha',0.2)
l1 = dat.cntrl{a}.latency(u);
l2 = dat.cool{a}.latency(u);
if ~isnan(l1)
    xline(l1,'k--')
end
if ~isnan(l2)
    xline(l2,'c--')
end
xlim([-1 2])
xlabel('time (sec)')
ylim([0 maxH])
if sp == 1
    ylabel('firing rate')
end
box on
sp = sp+1;
end

x = ages;
% x = agesJit;
% i = find(agIdx==3|agIdx==4);
% xA = [x(i);i];
% [~,sIdx] = sort(xA(1,:));
% xA = xA(:,sIdx);
% xA(1,:) = (0:5/(size(xA,2)-1):5)+ageGroups{3}(1);
% x(xA(2,:)) = xA(1,:);

minN = 20;
xL = [26 48];
xT = [25,30,35,40];
xTlbl = {'25','30','35','40+'};

subplot(3,2,3);hold on
sclSem = 100;
pVal = rPref.stats.sRank.p<0.01;
y1 = rPref.cntrl.ave;
sem1 = rPref.cntrl.sem;
n1 = rPref.cntrl.n;
scatter(x(n1>=minN&pVal),y1(n1>=minN&pVal),sem1(n1>=minN&pVal)*sclSem,'ko','MarkerFaceColor','flat','MarkerFaceAlpha',0.2,'LineWidth',1)
scatter(x(n1>=minN&~pVal),y1(n1>=minN&~pVal),sem1(n1>=minN&~pVal)*sclSem,'ko','MarkerFaceColor','none','MarkerFaceAlpha',0.2,'LineWidth',1)
y2 = rPref.cool.ave;
sem2 = rPref.cool.sem;
n2 = rPref.cool.n;
scatter(x(n2>=minN&pVal),y2(n2>=minN&pVal),sem2(n2>=minN&pVal)*sclSem,'co','MarkerFaceColor','flat','MarkerFaceAlpha',0.2,'LineWidth',1)
scatter(x(n2>=minN&~pVal),y2(n2>=minN&~pVal),sem2(n2>=minN&~pVal)*sclSem,'co','MarkerFaceColor','none','MarkerFaceAlpha',0.2,'LineWidth',1)
plot([x(n1>=minN);x(n2>=minN)],[y1(n1>=minN);y2(n2>=minN)],'k')
xlim(xL)
xticks(xT)
xticklabels(xTlbl)
ylabel('firing rate')
box on

subplot(3,2,4);hold on
sclSem = 5000;
pVal = si.stats.sRank.p<0.01;
y3 = si.ave;
sem3 = si.sem;
n3 = si.n;
yline(0,'k--')
scatter(x(n3>=minN&pVal),y3(n3>=minN&pVal),sem3(n3>=minN&pVal)*sclSem,'ko','MarkerFaceColor','flat','MarkerFaceAlpha',0.2,'LineWidth',1)
scatter(x(n3>=minN&~pVal),y3(n3>=minN&~pVal),sem3(n3>=minN&~pVal)*sclSem,'ko','MarkerFaceColor','none','MarkerFaceAlpha',0.2,'LineWidth',1)
ylim([-1,0.1])
xlim(xL)
xticks(xT)
xticklabels(xTlbl)
ylabel('SI')
box on

subplot(3,2,5);hold on
sclSem = 5000;
pVal = late.stats.sRank.p<0.01;
y4 = late.cntrl.ave;
sem4 = late.cntrl.sem;
n4 = late.cntrl.n;
scatter(x(n4>=minN&pVal),y4(n4>=minN&pVal),sem4(n4>=minN&pVal)*sclSem,'ko','MarkerFaceColor','flat','MarkerFaceAlpha',0.2,'LineWidth',1)
scatter(x(n4>=minN&~pVal),y4(n4>=minN&~pVal),sem4(n4>=minN&~pVal)*sclSem,'ko','MarkerFaceColor','none','MarkerFaceAlpha',0.2,'LineWidth',1)
y5 = late.cool.ave;
sem5 = late.cool.sem;
n5 = late.cool.n;
scatter(x(n5>=minN&pVal),y5(n5>=minN&pVal),sem5(n5>=minN&pVal)*sclSem,'co','MarkerFaceColor','flat','MarkerFaceAlpha',0.2,'LineWidth',1)
scatter(x(n5>=minN&~pVal),y5(n5>=minN&~pVal),sem5(n5>=minN&~pVal)*sclSem,'co','MarkerFaceColor','none','MarkerFaceAlpha',0.2,'LineWidth',1)
i = n4>=minN & n5>=minN;
plot([x(i);x(i)],[y4(i);y5(i)],'k')
xlim(xL)
xticks(xT)
xticklabels(xTlbl)
xlabel('age (postnatal day)')
ylabel('response latency (sec)')
box on

subplot(3,2,6);hold on
sclSem = 2000;
pVal = dLate.stats.sRank.p<0.01;
y6 = dLate.ave;
sem6 = dLate.sem;
n6 = dLate.n;
yline(0,'k--')
scatter(x(n6>=minN&pVal),y6(n6>=minN&pVal),sem6(n6>=minN&pVal)*sclSem,'ko','MarkerFaceColor','flat','MarkerFaceAlpha',0.2,'LineWidth',1)
scatter(x(n6>=minN&~pVal),y6(n6>=minN&~pVal),sem6(n6>=minN&~pVal)*sclSem,'ko','MarkerFaceColor','none','MarkerFaceAlpha',0.2,'LineWidth',1)
xlim(xL)
xticks(xT)
xticklabels(xTlbl)
xlabel('age (postnatal day)')
ylabel('delta response latency')
box on

%% plot groups

figure;
colororder({'k','c'})
D = vertcat(unitDataCntrl,unitDataCool);
idx = repmat(unitDataCntrl.good,2,1);
D = D(idx,:);
boxchart(D.AG,log10(D.rPref+0.1),'GroupByColor',D.manip,'Notch','on')
box on
% set(gca,'YScale','log')
% anovan(log10(D.rPref+0.1),{D.AG D.manip},'model','interaction','varnames',{'ageGroup','manip'})

figure;
D = vertcat(unitDataCntrl,unitDataCool);
idx = [unitDataCntrl.good;unitDataCool.good];
D = D(idx,:);
boxchart(D.AG,D.latency,'GroupByColor',D.manip)
% anovan(D.latency,{D.AG D.manip},'model','interaction','varnames',{'ageGroup','manip'})


%% split distribution

D = unitDataCntrl;
idx1 = unitDataCntrl.good & ~isnan(unitDataCntrl.latency); %pass during control
idx2 = unitDataCool.good & ~isnan(unitDataCool.latency); %pass during cooling
idx3 = idx1&idx2; %V1 independent
idx4 = idx1&(~idx2); %V1 dependent
idx5 = (~idx1)&idx2;

Y = D.latency;

figure; tiledlayout(length(ageGroups)+1,2)
nexttile; hold on
cdfplot(Y(idx1))
cdfplot(Y(idx3))
cdfplot(Y(idx4))
xticks([0:0.2:1])
xlim([0 1])
xlabel('latency (s)')
yticks([0 0.5 1])
ylabel('percentile')
title('All ages')
box on
axis square

X = vertcat([Y(idx3),ones(sum(idx3),1),zeros(sum(idx3),1),zeros(sum(idx3),1)],...
            [Y(idx4),ones(sum(idx4),1)+1,zeros(sum(idx4),1),ones(sum(idx4),1)]);

nexttile; hold on
bins = [0:0.05:1];
histogram(Y(idx1),bins)
histogram(Y(idx3),bins)
histogram(Y(idx4),bins)
xlabel('latency (s)')
ylabel('count')
legend({['pass cntrl; n=' num2str(sum(idx1))],...
        ['pass cntrl & cool; n=' num2str(sum(idx3))],...
        ['pass cntrl & ~cool; n=' num2str(sum(idx4))]},'Location','northeast')
box on
axis square

for ag = 1:length(ageGroups)
    nexttile; hold on
    cdfplot(Y(idx1&D.AG==ag))
    cdfplot(Y(idx3&D.AG==ag))
    cdfplot(Y(idx4&D.AG==ag))
    xticks([0:0.2:1])
    xlim([0 1])
    xlabel('latency (s)')
    yticks([0 0.5 1])
    ylabel('percentile')
    title(['Age group: P' num2str(ageGroups{ag}(1)) '-' num2str(ageGroups{ag}(2))])
    box on
    axis square
    
    nexttile; hold on
    histogram(Y(idx1&D.AG==ag),bins)
    histogram(Y(idx3&D.AG==ag),bins)
    histogram(Y(idx4&D.AG==ag),bins)
    xlabel('latency (s)')
    ylabel('count')
    legend({['pass cntrl; n=' num2str(sum(idx1&D.AG==ag))],...
            ['pass cntrl & cool; n=' num2str(sum(idx3&D.AG==ag))],...
            ['pass cntrl & ~cool; n=' num2str(sum(idx4&D.AG==ag))]},'Location','northeast')
    box on
    axis square

    Xtmp = vertcat([Y(idx3&D.AG==ag),repmat(((2*ag)-1)+2,sum(idx3&D.AG==ag),1),repmat(ag,sum(idx3&D.AG==ag),1),zeros(sum(idx3&D.AG==ag),1)],...
                   [Y(idx4&D.AG==ag),repmat((2*ag)+2,sum(idx4&D.AG==ag),1),repmat(ag,sum(idx4&D.AG==ag),1),ones(sum(idx4&D.AG==ag),1)]);
    X = vertcat(X,Xtmp);
end

sgtitle('control MU response latency')


% D = unitDataCool;
% idx1 = unitDataCool.good;
% idx2 = unitDataCntrl.good;
% idx3 = idx1&idx2;
% idx4 = idx1&(~idx2);
% idx5 = (~idx1)&idx2;
% 
% figure; tiledlayout(5,2)
% nexttile; hold on
% cdfplot(D.latency(idx1))
% cdfplot(D.latency(idx3))
% cdfplot(D.latency(idx4))
% xticks([0:0.2:1])
% xlim([0 1])
% xlabel('latency (s)')
% yticks([0 0.5 1])
% ylabel('percentile')
% title('All ages')
% 
% nexttile; hold on
% histogram(D.latency(idx1),[0:0.05:1])
% histogram(D.latency(idx3),[0:0.05:1])
% histogram(D.latency(idx4),[0:0.05:1])
% xlabel('latency (s)')
% ylabel('count')
% legend({['pass cool; n=' num2str(sum(idx1))],...
%         ['pass cool & cntrl; n=' num2str(sum(idx3))],...
%         ['pass cool & ~cntrl; n=' num2str(sum(idx4))]},'Location','northeast')
% 
% for ag = 1:4
%     nexttile; hold on
%     cdfplot(D.latency(idx1&D.AG==ag))
%     cdfplot(D.latency(idx3&D.AG==ag))
%     cdfplot(D.latency(idx4&D.AG==ag))
%     xticks([0:0.2:1])
%     xlim([0 1])
%     xlabel('latency (s)')
%     yticks([0 0.5 1])
%     ylabel('percentile')
%     title(['Age group:' num2str(ag)])
%     
%     nexttile; hold on
%     histogram(D.latency(idx1&D.AG==ag),[0:0.05:1])
%     histogram(D.latency(idx3&D.AG==ag),[0:0.05:1])
%     histogram(D.latency(idx4&D.AG==ag),[0:0.05:1])
%     xlabel('latency (s)')
%     ylabel('count')
%     legend({['pass cool; n=' num2str(sum(idx1&D.AG==ag))],...
%             ['pass cool & cntrl; n=' num2str(sum(idx3&D.AG==ag))],...
%             ['pass cool & ~cntrl; n=' num2str(sum(idx4&D.AG==ag))]},'Location','northeast')
% end
% 
% sgtitle('cool V1 MU response latency')
