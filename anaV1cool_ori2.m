%anaV1cool_ori2
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
% load(['/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/cooling/' proj '/V1cool_' anaMode '_ori_projectTbl.mat'])

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

%% Calculate Metrics

for a = 1:length(animals)
    if ~(isempty(dat.cntrl{a,1}) || isempty(dat.cool{a,1}))

        i = goodIdCntrl{a,1};

        r1 = dat.cntrl{a,1}.rPref(i);
        rPref.cntrl.dist{a} = r1;
        rPref.cntrl.ave(a) = mean(r1,'omitnan');
        rPref.cntrl.sem(a) = std(r1,'omitnan')/sqrt(length(r1));
        rPref.cntrl.n(a) = length(r1);

        r2 = dat.cool{a,1}.rPref(i);
        rPref.cool.dist{a} = r2;
        rPref.cool.ave(a) = mean(r2,'omitnan');
        rPref.cool.sem(a) = std(r2,'omitnan')/sqrt(length(r2));
        rPref.cool.n(a) = length(r2);

        si.dist{a} = (r2-r1)./(r2+r1);
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

        i = (goodIdCntrl{a,1}&goodIdCool{a,1}) & ~(isnan(dat.cntrl{a,1}.latency)|isnan(dat.cool{a,1}.latency));
        
        l1 = dat.cntrl{a,1}.latency(i);
        late.cntrl.dist{a} = l1;
        late.cntrl.ave(a) = mean(l1,'omitnan');
        late.cntrl.sem(a) = std(l1,'omitnan')/sqrt(length(l1));
        late.cntrl.n(a) = length(l1);

        l2 = dat.cool{a,1}.latency(i);
        late.cool.dist{a} = l2;
        late.cool.ave(a) = mean(l2,'omitnan');
        late.cool.sem(a) = std(l2,'omitnan')/sqrt(length(l2));
        late.cool.n(a) = length(l2);

        dLate.dist{a} = l2-l1;
        dLate.ave(a) = mean(dLate.dist{a},'omitnan');
        dLate.sem(a) = std(dLate.dist{a},'omitnan')/sqrt(length(dLate.dist{a}));
        dLate.n(a) = length(dLate.dist{a});

        if ~(isempty(l1)||isempty(l2))
            late.stats.sRank.p(a) = signrank(l1-l2);
            dLate.stats.sRank.p(a) = signrank(dLate.dist{a});
        else
            late.stats.sRank.p(a) = nan;
            dLate.stats.sRank.p(a) = nan;
        end

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

        late.cntrl.dist{a} = [];
        late.cntrl.ave(a) = nan;
        late.cntrl.sem(a) = nan;
        late.cntrl.n(a) = nan;

        late.cool.dist{a} = [];
        late.cool.ave(a) = nan;
        late.cool.sem(a) = nan;
        late.cool.n(a) = nan;

        dLate.dist{a} = [];
        dLate.ave(a) = nan;
        dLate.sem(a) = nan;
        dLate.n(a) = nan;

        late.stats.sRank.p(a) = nan;
        dLate.stats.sRank.p(a) = nan;

    end
end


%% Plot

figure('Position',[100 100 1000 700])

a = 4; 
count = 1;
for u = [42 52 54]
subplot(3,3,count); hold on
d1 = dat.cntrl{a}.spkTimes{u}(1,:);
d2 = dat.cool{a}.spkTimes{u}(1,:);
nTrials1 = max(dat.cntrl{a}.fr(u).trialNum,[],'all');
nTrials2 = max(dat.cool{a}.fr(u).trialNum,[],'all');
bw = 0.01;
bins = -1:bw:2;
h1 = histcounts(d1,bins)/(nTrials1*bw);
h2 = histcounts(d2,bins)/(nTrials2*bw);
maxH = 60;%max([h1 h2],[],'all')+1;
patch([0 1 1 0],[0 0 maxH maxH],'k','EdgeColor','none','FaceAlpha',0.2)
plot(bins(2:end),h1,'k')
plot(bins(2:end),h2,'c')
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
if count == 1
    ylabel('firing rate')
end
box on
count = count+1;
end

x = agesJit;
i = find(agIdx==3|agIdx==4);
xA = [x(i);i];
[~,sIdx] = sort(xA(1,:));
xA = xA(:,sIdx);
xA(1,:) = (0:5/(size(xA,2)-1):5)+ageGroups{3}(1);
x(xA(2,:)) = xA(1,:);

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


% subplot(3,2,3); hold on
% minN = 20;
% sclSem = 100;
% i = agIdx~=max(agIdx) & rPref.cntrl.n>=minN; 
% y1 = rPref.cntrl.ave(i); sem1 = rPref.cntrl.sem(i)*sclSem;
% y2 = rPref.cool.ave(i); sem2 = rPref.cool.sem(i)*sclSem;
% x = double(ages(i));
% for u = unique(x)
%     if sum(x==u)>1
%         x(x==u) = x(x==u)+(-0.3:0.6/(sum(x==u)-1):0.3);
%     end
% end
% scatter(x,y1,sem1,'ko','MarkerFaceColor','flat','MarkerFaceAlpha',0.2,'LineWidth',1)
% scatter(x,y2,sem2,'co','MarkerFaceColor','flat','MarkerFaceAlpha',0.2,'LineWidth',1)
% plot([x;x],[y1;y2],'k')
% i = agIdx~=max(agIdx) & rPref.cntrl.n<minN;
% y1 = rPref.cntrl.ave(i); sem1 = rPref.cntrl.sem(i)*sclSem;
% y2 = rPref.cool.ave(i);sem2 = rPref.cool.sem(i)*sclSem;
% x = double(ages(i));
% for u = unique(x)
%     if sum(x==u)>1
%         x(x==u) = x(x==u)+(-0.3:0.6/(sum(x==u)-1):0.3);
%     end
% end
% scatter(x,y1,sem1,'ko','MarkerFaceColor','none','MarkerFaceAlpha',0.2,'LineWidth',1)
% scatter(x,y2,sem2,'co','MarkerFaceColor','none','MarkerFaceAlpha',0.2,'LineWidth',1)
% plot([x;x],[y1;y2],'k')
% i = agIdx==max(agIdx) & rPref.cntrl.n>=minN;
% y1 = rPref.cntrl.ave(i); sem1 = rPref.cntrl.sem(i)*sclSem;
% y2 = rPref.cool.ave(i); sem2 = rPref.cool.sem(i)*sclSem;
% [~,sIdx] = sort(ages(i)); y1 = y1(sIdx); sem1 = sem1(sIdx); y2 = y2(sIdx); sem2 = sem2(sIdx);
% x = (0:5/length(y1):5)+41; x = x(2:end);
% scatter(x,y1,sem1,'ko','MarkerFaceColor','flat','MarkerFaceAlpha',0.2,'LineWidth',1)
% scatter(x,y2,sem2,'co','MarkerFaceColor','flat','MarkerFaceAlpha',0.2,'LineWidth',1)
% plot([x;x],[y1;y2],'k')
% i = agIdx==max(agIdx) & rPref.cntrl.n<minN;
% y1 = rPref.cntrl.ave(i); sem1 = rPref.cntrl.sem(i)*sclSem;
% y2 = rPref.cool.ave(i); sem2 = rPref.cool.sem(i)*sclSem;
% [~,sIdx] = sort(ages(i)); y1 = y1(sIdx); sem1 = sem1(sIdx); y2 = y2(sIdx); sem2 = sem2(sIdx);
% x = (0:5/length(y1):5)+41; x = x(2:end);
% scatter(x,y1,sem1,'ko','MarkerFaceColor','none','MarkerFaceAlpha',0.2,'LineWidth',1)
% scatter(x,y2,sem2,'co','MarkerFaceColor','none','MarkerFaceAlpha',0.2,'LineWidth',1)
% plot([x;x],[y1;y2],'k')
% xlim([25 50])
% ylabel('firing rate')
% box on
% 
% subplot(3,2,4); hold on
% minN = 20;
% sclSem = 5000;
% yline(0,'k--')
% i = agIdx~=max(agIdx) & si.n>=minN;
% y = si.ave(i); sem = si.sem(i)*sclSem;
% x = double(ages(i));
% for u = unique(x)
%     if sum(x==u)>1
%         x(x==u) = x(x==u)+(-0.3:0.6/(sum(x==u)-1):0.3);
%     end
% end
% scatter(x,y,sem,'ko','MarkerFaceColor','flat','MarkerFaceAlpha',0.2,'LineWidth',1)
% i = agIdx~=max(agIdx) & si.n<minN;
% y = si.ave(i); sem = si.sem(i)*sclSem;
% x = double(ages(i));
% for u = unique(x)
%     if sum(x==u)>1
%         x(x==u) = x(x==u)+(-0.3:0.6/(sum(x==u)-1):0.3);
%     end
% end
% scatter(x,y,sem,'ko','MarkerFaceColor','none','MarkerFaceAlpha',0.2,'LineWidth',1)
% i = agIdx==max(agIdx) & si.n>=minN;
% y = si.ave(i); sem = si.sem(i)*sclSem;
% [~,sIdx] = sort(ages(i)); y = y(sIdx); sem = sem(sIdx);
% x = (0:5/length(y):5)+41; x = x(2:end);
% scatter(x,y,sem,'ko','MarkerFaceColor','flat','MarkerFaceAlpha',0.2,'LineWidth',1)
% i = agIdx==max(agIdx) & si.n<minN;
% y = si.ave(i); sem = si.sem(i)*sclSem;
% [~,sIdx] = sort(ages(i)); y = y(sIdx); sem = sem(sIdx);
% x = (0:5/length(y):5)+41; x = x(2:end);
% scatter(x,y,sem,'ko','MarkerFaceColor','none','MarkerFaceAlpha',0.2,'LineWidth',1)
% xlim([25 50])
% ylim([-1 .1])
% ylabel('SI')
% box on
% 
% subplot(3,2,5); hold on
% minN = 20;
% sclSem = 5000;
% i = agIdx~=max(agIdx) & late.cntrl.n>=minN; 
% y1 = late.cntrl.ave(i); sem1 = late.cntrl.sem(i)*sclSem;
% y2 = late.cool.ave(i); sem2 = late.cool.sem(i)*sclSem;
% x = double(ages(i));
% for u = unique(x)
%     if sum(x==u)>1
%         x(x==u) = x(x==u)+(-0.3:0.6/(sum(x==u)-1):0.3);
%     end
% end
% scatter(x,y1,sem1,'ko','MarkerFaceColor','flat','MarkerFaceAlpha',0.2,'LineWidth',1)
% scatter(x,y2,sem2,'co','MarkerFaceColor','flat','MarkerFaceAlpha',0.2,'LineWidth',1)
% plot([x;x],[y1;y2],'k')
% i = agIdx~=max(agIdx) & late.cntrl.n<minN;
% y1 = late.cntrl.ave(i); sem1 = late.cntrl.sem(i)*sclSem;
% y2 = late.cool.ave(i);sem2 = late.cool.sem(i)*sclSem;
% x = double(ages(i));
% for u = unique(x)
%     if sum(x==u)>1
%         x(x==u) = x(x==u)+(-0.3:0.6/(sum(x==u)-1):0.3);
%     end
% end
% scatter(x,y1,sem1,'ko','MarkerFaceColor','none','MarkerFaceAlpha',0.2,'LineWidth',1)
% scatter(x,y2,sem2,'co','MarkerFaceColor','none','MarkerFaceAlpha',0.2,'LineWidth',1)
% plot([x;x],[y1;y2],'k')
% i = agIdx==max(agIdx) & late.cntrl.n>=minN;
% y1 = late.cntrl.ave(i); sem1 = late.cntrl.sem(i)*sclSem;
% y2 = late.cool.ave(i); sem2 = late.cool.sem(i)*sclSem;
% [~,sIdx] = sort(ages(i)); y1 = y1(sIdx); sem1 = sem1(sIdx); y2 = y2(sIdx); sem2 = sem2(sIdx);
% x = (0:5/length(y1):5)+41; x = x(2:end);
% scatter(x,y1,sem1,'ko','MarkerFaceColor','flat','MarkerFaceAlpha',0.2,'LineWidth',1)
% scatter(x,y2,sem2,'co','MarkerFaceColor','flat','MarkerFaceAlpha',0.2,'LineWidth',1)
% plot([x;x],[y1;y2],'k')
% i = agIdx==max(agIdx) & late.cntrl.n<minN;
% y1 = late.cntrl.ave(i); sem1 = late.cntrl.sem(i)*sclSem;
% y2 = late.cool.ave(i); sem2 = late.cool.sem(i)*sclSem;
% [~,sIdx] = sort(ages(i)); y1 = y1(sIdx); sem1 = sem1(sIdx); y2 = y2(sIdx); sem2 = sem2(sIdx);
% x = (0:5/length(y1):5)+41; x = x(2:end);
% scatter(x,y1,sem1,'ko','MarkerFaceColor','none','MarkerFaceAlpha',0.2,'LineWidth',1)
% scatter(x,y2,sem2,'co','MarkerFaceColor','none','MarkerFaceAlpha',0.2,'LineWidth',1)
% plot([x;x],[y1;y2],'k')
% xlim([25 50])
% xlabel('age (postnatal day)')
% ylabel('response latency (sec)')
% box on
% 
% subplot(3,2,6); hold on
% minN = 20;
% sclSem = 2000;
% yline(0,'k--')
% i = agIdx~=max(agIdx) & dLate.n>=minN;
% y = dLate.ave(i); sem = dLate.sem(i)*sclSem;
% x = double(ages(i));
% for u = unique(x)
%     if sum(x==u)>1
%         x(x==u) = x(x==u)+(-0.3:0.6/(sum(x==u)-1):0.3);
%     end
% end
% scatter(x,y,sem,'ko','MarkerFaceColor','flat','MarkerFaceAlpha',0.2,'LineWidth',1)
% i = agIdx~=max(agIdx) & dLate.n<minN;
% y = dLate.ave(i); sem = dLate.sem(i)*sclSem;
% x = double(ages(i));
% for u = unique(x)
%     if sum(x==u)>1
%         x(x==u) = x(x==u)+(-0.3:0.6/(sum(x==u)-1):0.3);
%     end
% end
% scatter(x,y,sem,'ko','MarkerFaceColor','none','MarkerFaceAlpha',0.2,'LineWidth',1)
% i = agIdx==max(agIdx) & dLate.n>=minN;
% y = dLate.ave(i); sem = dLate.sem(i)*sclSem;
% [~,sIdx] = sort(ages(i)); y = y(sIdx); sem = sem(sIdx);
% x = (0:5/length(y):5)+41; x = x(2:end);
% scatter(x,y,sem,'ko','MarkerFaceColor','flat','MarkerFaceAlpha',0.2,'LineWidth',1)
% i = agIdx==max(agIdx) & dLate.n<minN;
% y = dLate.ave(i); sem = dLate.sem(i)*sclSem;
% [~,sIdx] = sort(ages(i)); y = y(sIdx); sem = sem(sIdx);
% x = (0:5/length(y):5)+41; x = x(2:end);
% scatter(x,y,sem,'ko','MarkerFaceColor','none','MarkerFaceAlpha',0.2,'LineWidth',1)
% xlim([25 50])
% xlabel('age (postnatal day)')
% ylabel('delta response latency')
% box on





% figure
% subplot(2,2,1); hold on
% idx = agIdx==1|agIdx==2;
% xLims = [min(ages(idx))-2 max(ages(idx))+2 max(ages(idx))+2 min(ages(idx))-2];
% ag3 = rPref.cntrl.ave(agIdx==3);
% mAg3 = mean(ag3,'omitnan');
% semAg3 = std(ag3,'omitnan')/sqrt(length(ag3));
% y1 = mAg3-semAg3;
% y2 = mAg3+semAg3;
% y = [y1 y1 y2 y2];
% patch(xLims,y,'k','EdgeColor','none','FaceAlpha',0.2)
% yline(mAg3,'k--')
% ag3 = rPref.cool.ave(agIdx==3);
% mAg3 = mean(ag3,'omitnan');
% semAg3 = std(ag3,'omitnan')/sqrt(length(ag3));
% y1 = mAg3-semAg3;
% y2 = mAg3+semAg3;
% y = [y1 y1 y2 y2];
% patch(xLims,y,'c','EdgeColor','none','FaceAlpha',0.2)
% yline(mAg3,'c--')
% errorbar(ages(idx),rPref.cntrl.ave(idx),rPref.cntrl.sem(idx),'k','LineStyle','none')
% errorbar(ages(idx),rPref.cool.ave(idx),rPref.cool.sem(idx),'c','LineStyle','none')
% xlabel('age')
% xlim(xLims(1:2))
% ylabel('firing rate')
% box on
% 
% subplot(2,2,2); hold on
% ag3 = si.ave(agIdx==3);
% mAg3 = mean(ag3,'omitnan');
% semAg3 = std(ag3,'omitnan')/sqrt(length(ag3));
% y1 = mAg3-semAg3;
% y2 = mAg3+semAg3;
% y = [y1 y1 y2 y2];
% patch(xLims,y,'r','EdgeColor','none','FaceAlpha',0.2)
% yline(mAg3,'r--')
% yline(0,'k--')
% errorbar(ages(idx),si.ave(idx),si.sem(idx),'LineStyle','none')
% xlabel('age')
% xlim(xLims(1:2))
% ylabel('si')
% box on
% 
% subplot(2,2,3); hold on
% ag3 = late.cntrl.ave(agIdx==3);
% mAg3 = mean(ag3,'omitnan');
% semAg3 = std(ag3,'omitnan')/sqrt(length(ag3));
% y1 = mAg3-semAg3;
% y2 = mAg3+semAg3;
% y = [y1 y1 y2 y2];
% patch(xLims,y,'k','EdgeColor','none','FaceAlpha',0.2)
% yline(mAg3,'k--')
% ag3 = late.cool.ave(agIdx==3);
% mAg3 = mean(ag3,'omitnan');
% semAg3 = std(ag3,'omitnan')/sqrt(length(ag3));
% y1 = mAg3-semAg3;
% y2 = mAg3+semAg3;
% y = [y1 y1 y2 y2];
% patch(xLims,y,'c','EdgeColor','none','FaceAlpha',0.2)
% yline(mAg3,'c--')
% errorbar(ages(idx),late.cntrl.ave(idx),late.cntrl.sem(idx),'k','LineStyle','none')
% errorbar(ages(idx),late.cool.ave(idx),late.cool.sem(idx),'c','LineStyle','none')
% xlabel('age')
% xlim(xLims(1:2))
% ylabel('latency')
% box on
% 
% subplot(2,2,4); hold on
% ag3 = dLate.ave(agIdx==3);
% mAg3 = mean(ag3,'omitnan');
% semAg3 = std(ag3,'omitnan')/sqrt(length(ag3));
% y1 = mAg3-semAg3;
% y2 = mAg3+semAg3;
% y = [y1 y1 y2 y2];
% patch(xLims,y,'r','EdgeColor','none','FaceAlpha',0.2)
% yline(mAg3,'r--')
% yline(0,'k--')
% errorbar(ages(idx),dLate.ave(idx),dLate.sem(idx),'LineStyle','none')
% xlabel('age')
% xlim(xLims(1:2))
% ylabel('delta latency')
% box on



%% Make data table

animalData1 = table();
animalData1.id = repmat(animals,2,1);
% animalData1.age = repmat(ages',2,1);
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
% animalData2.age = ages';
animalData2.ageJit = agesJit';
animalData2.SI_mean = si.ave';
animalData2.SI_sem = si.sem';
animalData2.SI_n = si.n';
animalData2.dL_mean = dLate.ave';
animalData2.dL_sem = dLate.sem';
animalData2.dL_n = dLate.n';

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
animalData.Lcntrl_mean = late.cntrl.ave';
animalData.Lcntrl_sem = late.cntrl.sem';
animalData.Lcntrl_n = late.cntrl.n';
% animalData.Lcntrl_dist = late.cntrl.dist';
animalData.Lcool_mean = late.cool.ave';
animalData.Lcool_sem = late.cool.sem';
animalData.Lcool_n = late.cool.n';
% animalData.Lcool_dist = late.cool.dist';
animalData.dL_mean = dLate.ave';
animalData.dL_sem = dLate.sem';
animalData.dL_n = dLate.n';
% animalData.dL_dist = dLate.dist';

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
unitDataCool.good = screenUnits(cool,'MU');
clear id uAge uAG manip g

unitData = vertcat(unitDataCntrl,unitDataCool);
unitData = unitData(unitData.good,:);

%% split distribution

D = unitDataCntrl;
idx1 = unitDataCntrl.good;
idx2 = unitDataCool.good;
idx3 = idx1&idx2;
idx4 = idx1&(~idx2);
idx5 = (~idx1)&idx2;

figure; tiledlayout(5,2)
nexttile; hold on
cdfplot(D.latency(idx1))
cdfplot(D.latency(idx3))
cdfplot(D.latency(idx4))
xticks([0:0.2:1])
xlim([0 1])
xlabel('latency (s)')
yticks([0 0.5 1])
ylabel('percentile')
title('All ages')

nexttile; hold on
histogram(D.latency(idx1),[0:0.05:1])
histogram(D.latency(idx3),[0:0.05:1])
histogram(D.latency(idx4),[0:0.05:1])
xlabel('latency (s)')
ylabel('count')
legend({['pass cntrl; n=' num2str(sum(idx1))],...
        ['pass cntrl & cool; n=' num2str(sum(idx3))],...
        ['pass cntrl & ~cool; n=' num2str(sum(idx4))]},'Location','northeast')

for ag = 1:4
    nexttile; hold on
    cdfplot(D.latency(idx1&D.AG==ag))
    cdfplot(D.latency(idx3&D.AG==ag))
    cdfplot(D.latency(idx4&D.AG==ag))
    xticks([0:0.2:1])
    xlim([0 1])
    xlabel('latency (s)')
    yticks([0 0.5 1])
    ylabel('percentile')
    title(['Age group:' num2str(ag)])
    
    nexttile; hold on
    histogram(D.latency(idx1&D.AG==ag),[0:0.05:1])
    histogram(D.latency(idx3&D.AG==ag),[0:0.05:1])
    histogram(D.latency(idx4&D.AG==ag),[0:0.05:1])
    xlabel('latency (s)')
    ylabel('count')
    legend({['pass cntrl; n=' num2str(sum(idx1&D.AG==ag))],...
            ['pass cntrl & cool; n=' num2str(sum(idx3&D.AG==ag))],...
            ['pass cntrl & ~cool; n=' num2str(sum(idx4&D.AG==ag))]},'Location','northeast')
end

sgtitle('control MU response latency')


D = unitDataCool;
idx1 = unitDataCool.good;
idx2 = unitDataCntrl.good;
idx3 = idx1&idx2;
idx4 = idx1&(~idx2);
idx5 = (~idx1)&idx2;

figure; tiledlayout(5,2)
nexttile; hold on
cdfplot(D.latency(idx1))
cdfplot(D.latency(idx3))
cdfplot(D.latency(idx4))
xticks([0:0.2:1])
xlim([0 1])
xlabel('latency (s)')
yticks([0 0.5 1])
ylabel('percentile')
title('All ages')

nexttile; hold on
histogram(D.latency(idx1),[0:0.05:1])
histogram(D.latency(idx3),[0:0.05:1])
histogram(D.latency(idx4),[0:0.05:1])
xlabel('latency (s)')
ylabel('count')
legend({['pass cool; n=' num2str(sum(idx1))],...
        ['pass cool & cntrl; n=' num2str(sum(idx3))],...
        ['pass cool & ~cntrl; n=' num2str(sum(idx4))]},'Location','northeast')

for ag = 1:4
    nexttile; hold on
    cdfplot(D.latency(idx1&D.AG==ag))
    cdfplot(D.latency(idx3&D.AG==ag))
    cdfplot(D.latency(idx4&D.AG==ag))
    xticks([0:0.2:1])
    xlim([0 1])
    xlabel('latency (s)')
    yticks([0 0.5 1])
    ylabel('percentile')
    title(['Age group:' num2str(ag)])
    
    nexttile; hold on
    histogram(D.latency(idx1&D.AG==ag),[0:0.05:1])
    histogram(D.latency(idx3&D.AG==ag),[0:0.05:1])
    histogram(D.latency(idx4&D.AG==ag),[0:0.05:1])
    xlabel('latency (s)')
    ylabel('count')
    legend({['pass cool; n=' num2str(sum(idx1&D.AG==ag))],...
            ['pass cool & cntrl; n=' num2str(sum(idx3&D.AG==ag))],...
            ['pass cool & ~cntrl; n=' num2str(sum(idx4&D.AG==ag))]},'Location','northeast')
end

sgtitle('cool V1 MU response latency')
