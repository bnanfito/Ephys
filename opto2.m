clear all
close all

anaMode = 'MU';
dataFold = 'F:\Brandon\data';
animal = 'febs6'; unit = '000'; expt = '003'; probe = 1;
baseName = [animal '_u' unit '_' expt];

physDir = fullfile(dataFold,'Ephys',animal,baseName);
load(fullfile(physDir,[baseName '_trialInfo.mat']))
load(fullfile(physDir,[baseName '.analyzer']),'-mat','Analyzer')
load(fullfile(physDir,[baseName '_id.mat']))

predelay = getparam('predelay',Analyzer);
stimTime = getparam('stim_time',Analyzer);
postdelay = getparam('postdelay',Analyzer);
stimStartID = 3;
tLims = [1,stimTime];

switch anaMode
    case 'SU'
SUTrialData(fullfile(dataFold,'Ephys'),animal,unit,expt,probe,'id',stimStartID,tLims(1),tLims(2),0)
load(fullfile(physDir,[baseName '_p' num2str(probe) '_SUTrial.mat']),'SU','SUinfo')
uDat = SU; uInfo = SUinfo;
clear SU SUinfo
    case 'MU'
MUThreshTrialData(fullfile(dataFold,'Ephys'),animal,unit,expt,probe,'id',stimStartID,tLims(1),tLims(2),0)
load(fullfile(physDir,[baseName '_p' num2str(probe) '_MUThreshTrial.mat']),'MUThresh','MUThreshInfo')
uDat = MUThresh; uInfo = MUThreshInfo;
clear MUThresh MUThreshInfo
end

nU = size(uDat,2);
nConds = size(trialInfo.domval,1);
nTrials = length(trialInfo.triallist);

clrs = {'k','r','b'};

oriIdx = strcmp(trialInfo.dom,'ori');
oris = unique(trialInfo.domval(:,oriIdx));
oriDom = trialInfo.domval(:,oriIdx);
optoIdx = strcmp(trialInfo.dom,'useFuncGen');
optoDom = trialInfo.domval(:,optoIdx);
dcIdx = strcmp(trialInfo.dom,'dutyCycle_C1');
dcDom = trialInfo.domval(:,dcIdx);
if ~isempty(trialInfo.blankId)
    oriDom = [oriDom;nan];
    optoDom = [optoDom;nan];
    dcDom = [dcDom;nan];
end
trialOris = oriDom(trialInfo.triallist);
trialOpto = optoDom(trialInfo.triallist);
trialDC = dcDom(trialInfo.triallist);
for u = 1:nU
    figure
    tl = tiledlayout('flow','TileSpacing','compact');
    bcfr = uDat(u).stimFrate-mean(uDat(u).baseFrate);
    for o = 1:length(oris)
        rCntrl{u}(:,o) = bcfr(trialOris==oris(o)&trialOpto==0);
        rOptoConst{u}(:,o) = bcfr(trialOris==oris(o)&trialOpto==1&trialDC==100);
        rOptoPulse{u}(:,o) = bcfr(trialOris==oris(o)&trialOpto==1&trialDC==50);
    end

    rasterX = [];
    rasterY = [];
    for t = 1:nTrials
        curSpks = uDat(u).spktimes{t}/1000;
        rasterX = [rasterX curSpks];
        rasterY = [rasterY repmat(t,1,length(curSpks))];
        [sdfY{u}(t,:),sdfX] = compute_sdf2(curSpks,0.001,[-tLims(1) tLims(2)]);
    end

    rMaxCntrl(u) = max(mean(rCntrl{u}));
    rMaxOptoConst(u) = max(mean(rOptoConst{u}));
    rMaxOptoPulse(u) = max(mean(rOptoPulse{u}));

    ax1 = nexttile;hold on
    for c = 1:3
        switch c
            case 1
                curTrials = find(trialOpto==0);
            case 2
                curTrials = find(trialOpto==1 & trialDC==100);
            case 3
                curTrials = find(trialOpto==1 & trialDC==50);
        end
        lbl{c} = [trialInfo.dom{1} ' ' num2str(trialInfo.domval(c))];
        idx = ismember(rasterY,curTrials);
        plot(rasterX(idx),rasterY(idx),'.','Color',clrs{c})
    end
    
    nexttile
    x = deg2rad(oris); x = [x;x(1)];
    y = mean(rCntrl{u}); y = [y,y(1)];
    polarplot(x,y,clrs{1},'LineWidth',2);hold on
    y = mean(rOptoConst{u}); y = [y,y(1)];
    polarplot(x,y,clrs{2},'LineWidth',2)
    y = mean(rOptoPulse{u}); y = [y,y(1)];
    polarplot(x,y,clrs{3},'LineWidth',2)

    ax2 = nexttile;hold on
    for c = 1:3
        switch c
            case 1
                curTrials = find(trialOpto==0);
            case 2
                curTrials = find(trialOpto==1 & trialDC==100);
            case 3
                curTrials = find(trialOpto==1 & trialDC==50);
        end
        p(c) = plot(sdfX,mean(sdfY{u}(curTrials,:),'omitnan'),'LineWidth',2,'Color',clrs{c});
    end
    linkaxes([ax1 ax2],'x')
    lgd = legend(p,{'cntrl','opto const','opto pulse'});
    lgd.Layout.Tile = 4;
    title([baseName ' ' anaMode num2str(u)])

    figName = [baseName '_' anaMode num2str(u)];
%     saveas(gcf,fullfile(physDir,figName),'fig')
%     saveas(gcf,fullfile(physDir,figName),'svg')
%     saveas(gcf,fullfile(physDir,figName),'tif')
end

figure;scatter(rMaxCntrl,rMaxOptoConst,'r')
hold on;scatter(rMaxCntrl,rMaxOptoPulse,'b')
plot([-10 60],[-10 60],'k--')
xline(0,'k');yline(0,'k')
box on; axis square

figure;hold on
a = rMaxCntrl;b = rMaxOptoConst;
SIconst = (b-a)./(b+a);
h = histogram(SIconst,'BinWidth',0.1);
h.EdgeColor = 'none';
h.FaceColor = 'r';
h.FaceAlpha = 0.5;

a = rMaxCntrl;b = rMaxOptoPulse;
SIpulse = (b-a)./(b+a);
h = histogram(SIpulse,'BinWidth',0.1);
h.EdgeColor = 'none';
h.FaceColor = 'b';
h.FaceAlpha = 0.5;

xlabel('SI')
ylabel('count')



