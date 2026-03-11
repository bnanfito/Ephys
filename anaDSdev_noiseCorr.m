%anaDSdev_noiseCorr
clear all
close all

anaMode = 'MU';
load(['Y:\Brandon\data\dataSets\DSdev\DSdev_' anaMode 'dataSet.mat'])

sC = [];
nC = [];
pAge = [];
areas = {'V1','PSS'};
pArea = [];
eId = [];
for e = 1:height(projTbl)
    data = projTbl.sumStats{e};
    data = data(screenUnits(data,anaMode),:);
    if isempty(data) || height(data)==1
        continue
    end
    
    for u = 1:height(data)
        r{e,u} = data.response{u};
        rM{e}(u,:) = mean(r{e,u},'omitnan');
        rN{e,u} = r{e,u}-rM{e}(u,:);
    end
    
    sigCorr = nan(height(data));
    noiseCorr = nan(height(data));
    for uA = 1:height(data)
        for uB = uA:height(data)
            rN_A = rN{e,uA}(:);
            rN_B = rN{e,uB}(:);
            noiseCorr(uA,uB) = corr(rN_A,rN_B);
    
            rS_A = rM{e}(uA,:);
            rS_B = rM{e}(uB,:);
            sigCorr(uA,uB) = corr(rS_A',rS_B');
        end
    end
    sigCorr(eye(height(data))==1) = nan;
    sigCorr = sigCorr(~isnan(sigCorr));
    noiseCorr(eye(height(data))==1) = nan;
    noiseCorr = noiseCorr(~isnan(noiseCorr));

    pAge = [pAge ;repmat(projTbl.age(e),length(noiseCorr),1)];
    pArea = [pArea ;repmat( find(strcmp(projTbl.recSite{e},areas)) ,length(noiseCorr),1)];
    nC = [nC; noiseCorr];
    sC = [sC; sigCorr];
    eId = [eId; repmat(e,length(noiseCorr),1)];

    lm = fitlm(sigCorr,noiseCorr);
    slope(e) = lm.Coefficients.Estimate(2);
    slope_pval(e) = lm.Coefficients.pValue(2);
end

figure
ageGroups = {[29 32],[33 36],[37 52]};
for ar = 1:2
subplot(1,2,ar);hold on
for ag = 1:3
idx = pArea==ar & pAge>=ageGroups{ag}(1) & pAge<=ageGroups{ag}(2);
x = sC(idx);
y = nC(idx);
plot(x,y,'.')
lm = fitlm(x,y);
plot(x,lm.Fitted)
xlim([-1 1])
ylim([-1 1])
axis square
end
end


figure;hold on
idx = strcmp(projTbl.recSite,'V1')
plot(projTbl.age(idx),slope(idx),'o')
idx = strcmp(projTbl.recSite,'PSS')
plot(projTbl.age(idx),slope(idx),'o')
